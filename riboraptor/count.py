"""Utilities for read counting operations."""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

from collections import Counter
from collections import OrderedDict
import os
import pickle
import re
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
import pybedtools
from pyfaidx import Fasta
import pysam
from scipy.stats import norm

from .wig import WigReader
from .helpers import mkdir_p

__SAM_NOT_UNIQ_FLAGS__ = [4, 20, 256, 272, 2048]

# Unmapped, Unmapped+Reverse strand, Not primary alignment,
# Not primary alignment + reverse strand, supplementary alignment

# Source: https://broadinstitute.github.io/picard/explain-flags.html


def _check_file_exists(filepath):
    if os.path.isfile(os.path.abspath(filepath)):
        return True
    return False


def _create_bam_index(bam):
    """Create bam index.

    Parameters
    ----------
    bam : str
          Path to bam file
    """
    if not os.path.exists('{}.bai'.format(bam)):
        pysam.index(bam)


def _is_read_uniq_mapping(read):
    """Check if read is uniquely mappable.

    Parameters
    ----------

    read : pysam.Alignment.fetch object


    Most reliable: ['NH'] tag
    """
    # Filter out secondary alignments
    if read.is_secondary:
        return False
    tags = dict(read.get_tags())
    try:
        nh_count = tags['NH']
    except KeyError:
        # Reliable in case of STAR
        if read.mapping_quality == 255:
            return True
        if read.mapping_quality < 1:
            return False
        # NH tag not set so rely on flags
        if read.flag in __SAM_NOT_UNIQ_FLAGS__:
            return False
        else:
            raise RuntimeError('Malformed BAM?')
    if nh_count == 1:
        return True
    return False


def bedgraph_to_bigwig(bedgraph, chrom_sizes,
                       saveto, input_is_stream=False):
    """Convert bedgraph to bigwig.

    Parameters
    ----------
    bedgraph : str
               Path to bedgraph file
    chrom_sizes : str
                  Path to genome chromosome sizes file
    saveto : str
             Path to write bigwig file
    input_is_stream : bool
                      True if input is sent through stdin
    """
    if input_is_stream:
        with tempfile.TemporaryFile() as fp:
            fp.write(('\n').join(bedgraph))
            filename = fp.name
        bedgraph = filename

    cmds = ['bedSort', bedgraph, bedgraph]
    p = subprocess.Popen(cmds, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, universal_newlines=True)
    stdout, stderr = p.communicate()

    cmds = ['bedGraphToBigWig', bedgraph, chrom_sizes, saveto]
    p = subprocess.Popen(cmds, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, universal_newlines=True)
    stdout, stderr = p.communicate()


def bam_to_bedgraph(bam, strand='both', end_type='5prime', saveto=None):
    """Create bigwig from bam.

    Parameters
    ----------
    bam : str
          Path to bam file
    strand : str, optional
             Use reads mapping to '+/-/both' strands
    end_type : str
               Use only end_type=5prime(5') or "3prime(3')"
    saveto : str, optional
              Path to write bedgraph

    Returns
    -------
    genome_cov : str
                 Bedgraph output

    """
    if strand not in ['+', '-', 'both']:
        raise RuntimeError('Strand should be one of \'+\', \'-\', \'both\'')
    if end_type == '5prime':
        extra_args = '-5'
    elif end_type == '3prime':
        extra_args = '-3'
    elif end_type == 'either':
        extra_args = ''
    bed = pybedtools.BedTool(bam)
    if strand != 'both':
        genome_cov = bed.genome_coverage(bg=True,
                                         strand=strand,
                                         additional_args=extra_args)
    else:
        genome_cov = bed.genome_coverage(bg=True,
                                         additional_args=extra_args)
    print(type(genome_cov))
    if saveto:
        with open(saveto, 'w') as outf:
            outf.write(str(genome_cov))
    return genome_cov


def read_enrichment(read_lengths,
                    enrichment_range=range(28, 33),
                    input_is_stream=False,
                    input_is_file=False):
    """Calculate read enrichment for a certain range of lengths

    Parameters
    ----------
    read_lengths : Counter
                       A counter with read lengths and their counts
    enrichment_range : range or str
                       Range of reads to concentrate upon
                       (28-32 or range(28,33))
    input_is_stream : bool
                      True if input is sent through stdin

    Returns
    -------
    ratio : float
             Enrichment in this range

    """
    if input_is_file:
        if not _check_file_exists(read_lengths):
            raise RuntimeError('{} does not exist.'.format(read_lengths))
        read_lengths = pd.read_table(read_lengths,
                                     names=['frag_len', 'frag_count'],
                                     sep='\t')
        read_lengths = pd.Series(read_lengths.frag_count.tolist(), index=read_lengths.frag_len.tolist())
    elif input_is_stream:
        counter = {}
        for line in read_lengths:
            splitted = list(map(lambda x: int(x), line.strip().split('\t')))
            counter[splitted[0]] = splitted[1]
        read_lengths = Counter(counter)
    if isinstance(read_lengths, Counter):
        read_lengths = pd.Series(read_lengths)
    if isinstance(enrichment_range, unicode) or\
            isinstance(enrichment_range, str):
        splitted = list(
            map(lambda x: int(x), enrichment_range.strip().split('-')))
        enrichment_range = range(splitted[0], splitted[1]+1)
    rpf_signal = read_lengths[enrichment_range].sum()
    total_signal = read_lengths.sum()
    array = [[x] * y for x, y in sorted(read_lengths.iteritems())]
    mean_length, std_dev_length = norm.fit(
        np.concatenate(array).ravel().tolist())

    # mean_length_floor = np.floor(mean_length)
    # 1 - P(x1 < X <x2) = P(X<x1) + P(X>x2) = cdf(x1) + sf(x2)
    cdf_min = norm.cdf(min(enrichment_range), mean_length, std_dev_length)
    sf_max = norm.sf(max(enrichment_range), mean_length, std_dev_length)
    pvalue = cdf_min + sf_max
    ratio = rpf_signal / float(total_signal - rpf_signal)
    return ratio, pvalue


def gene_coverage(gene_name, bed, bw, master_offset=0):
    """Get gene coverage.

    Parameters
    ----------
    gene_name : str
                Gene name
    bed : str
          Path to CDS or 5'UTR or 3'UTR bed
    bw : str
         Path to bigwig to fetch the scores from
    master_offset : int
                    Number of bases to count upstream

    Returns
    -------
    coverage_combined : series
                        Series with index as position and value as coverage
    intervals_for_fasta_read : list
                                List of tuples
    index_to_genomic_pos_map : series

    gene_offset : int
                Gene wise offsets
    """
    bw = WigReader(bw)
    chromsome_lengths = bw.get_chromosomes
    bed = pybedtools.BedTool(bed).to_dataframe()
    assert gene_name in bed['name'].tolist()
    gene_group = bed[bed['name'] == gene_name]

    assert len(gene_group['strand'].unique()) == 1
    assert len(gene_group['chrom'].unique()) == 1
    chrom = gene_group['chrom'].unique()[0]
    strand = gene_group['strand'].unique()[0]

    # Collect all intervals at once
    intervals = zip(gene_group['chrom'], gene_group['start'],
                    gene_group['end'], gene_group['strand'])
    chrom_length = chromsome_lengths[str(intervals[0][0])]
    # Need to convert to list instead frm tuples
    # TODO fix this?
    intervals = map(list, intervals)
    if strand == '+':
        # For positive strand shift
        # start codon position intervals[0][1] by -offset
        if intervals[0][1] - master_offset >= 0:
            intervals[0][1] = intervals[0][1] - master_offset
            gene_offset = master_offset
        else:
            sys.stderr.write(
                'Cannot offset beyond 0 for interval: {}. \
                Set to start of chromsome.\n'.format(intervals[0]))
            # Reset offset to minimum possible
            gene_offset = intervals[0][1]
            intervals[0][1] = 0
    else:
        # Else shift cooridnate of last element in intervals stop by + offset
        if (intervals[-1][2] + master_offset <= chrom_length):
            intervals[-1][2] = intervals[-1][2] + master_offset
            gene_offset = master_offset
        else:
            sys.stderr.write('Cannot offset beyond 0 for interval: {}. \
                             Set to end of chromsome.\n'.format(intervals[-1]))
            gene_offset = chrom_length - intervals[-1][2]
            # 1-end so chrom_length
            intervals[-1][2] = chrom_length

    intervals = map(tuple, intervals)

    interval_coverage_list = []
    for index, coverage in enumerate(bw.read(intervals)):
        strand = intervals[index][3]
        if strand == '+':
            series_range = range(intervals[index][1], intervals[index][2])
        elif strand == '-':
            series_range = range(intervals[index][2], intervals[index][1], -1)

        series = pd.Series(coverage, index=series_range)
        interval_coverage_list.append(series)

    if len(interval_coverage_list) == 0:
        # Some genes might not be present in the bigwig at all
        sys.stderr.write(
            'Got empty list! intervals  for chr : {}\n'.format(intervals[0][0]))
        return ([], None)

    coverage_combined = interval_coverage_list[0]
    for interval_coverage in interval_coverage_list[1:]:
        coverage_combined = coverage_combined.combine_first(interval_coverage)
    coverage_combined = coverage_combined.fillna(0)
    coverage_index = np.arange(len(coverage_combined)) - gene_offset
    index_to_genomic_pos_map = pd.Series(coverage_combined.index.tolist(),
                                         index=coverage_index)
    intervals_for_fasta_read = []
    for pos in index_to_genomic_pos_map.values:
        intervals_for_fasta_read.append((chrom, pos, pos + 1, strand))
    coverage_combined = coverage_combined.reset_index(drop=True)
    coverage_combined = coverage_combined.rename(lambda x: x - gene_offset)

    return (coverage_combined, intervals_for_fasta_read,
            index_to_genomic_pos_map, gene_offset)


def get_fasta_sequence(fasta, intervals):
    """Extract fasta sequence given a list of intervals.

    Parameters
    ----------
    fasta : str
            Path to fasta file

    intervals : list(tuple)
                A list of tuple in the form [(chrom, start, stop, strand)]

    Returns
    -------
    seq : list
          List of sequences at intervals
    """
    fasta = Fasta(fasta)
    sequence = []
    for interval in intervals:
        chrom, start, stop, strand = interval
        if strand == '+':
            seq = fasta[chrom][int(start):int(stop)].seq
        elif strand == '-':
            seq = fasta[chrom][int(start):int(stop)].reverse.seq
        sequence.append(seq)
    return sequence


def get_region_sizes(bed):
    """Get collapsed lengths of gene in bed.

    Parameters
    ----------
    bed : str
          Path to bed file

    Returns
    -------
    region_sizes : dict
                   Region sies with gene names as key
                   and value as size of this named region
    """
    bed = pybedtools.BedTool(bed).to_dataframe()
    region_bed_grouped = bed.groupby('name')
    region_sizes = {}
    for gene_name, gene_group in region_bed_grouped:
        # Get rid of trailing dots
        # TODO: this can be letha for genomes with names
        # such as orf19.123
        # gene_name = re.sub(r'\.[0-9]+', '', gene_name)
        # Collect all intervals at once
        intervals = zip(gene_group['chrom'], gene_group['start'],
                        gene_group['end'], gene_group['strand'])
        for interval in intervals:
            if gene_name not in region_sizes:
                # End is always 1-based so does not require +1
                region_sizes[gene_name] = interval[2] - interval[1]
            else:
                region_sizes[gene_name] += interval[2] - interval[1]
    sizes_df = pd.DataFrame(region_sizes.items(), columns=['name', 'length'])
    sizes_df.sort_values(by='length', ascending=False, inplace=True)
    region_sizes = OrderedDict([tuple(x)
                                for x in sizes_df[['name', 'length']].values])
    return region_sizes


def htseq_to_cpm(htseq_f, saveto=None):
    """Convert HTSeq counts to CPM.

    Parameters
    ----------
    htseq_f : str
              Path to HTseq counts file
    saveto : str, optional
              Path to output file

    Returns
    -------
    cpm : dataframe
          CPM

    """
    htseq = read_htseq(htseq_f)
    rate = htseq['counts']
    denom = rate.sum()
    cpm = rate / denom * 1e6
    if saveto:
        pd.DataFrame(cpm, columns=['cpm']).to_csv(
            saveto, sep='\t', index=True, header=False)
    cpm = pd.DataFrame(cpm)  # , columns=['cpm'])
    cpm.columns = ['cpm']
    cpm = cpm.sort_values(by='cpm', ascending=False)
    return cpm


def htseq_to_tpm(htseq_f, cds_bed_f, saveto=None):
    """Convert HTSeq counts to TPM.

    Parameters
    ----------
    htseq_f : str
              Path to HTseq counts file
    region_sizes : dict
                   Dict with keys as gene and values
                   as length (CDS/Exon) of that gene
    saveto : str, optional
              Path to output file

    Returns
    -------
    tpm : dataframe
          TPM

    """
    cds_bed_sizes = get_region_sizes(cds_bed_f)
    htseq = read_htseq(htseq_f)
    rate = np.log(htseq['counts']).subtract(np.log(cds_bed_sizes))
    denom = np.log(np.sum(np.exp(rate)))
    tpm = np.exp(rate - denom + np.log(1e6))
    if saveto:
        pd.DataFrame(tpm, columns=['tpm']).to_csv(
            saveto, sep='\t', index=True, header=False)
    tpm = pd.DataFrame(tpm, columns=['tpm'])
    tpm = tpm.sort_values(by='tpm', ascending=False)
    return tpm


def mapping_reads_summary(bam):
    """Count number of mapped reads.

    Parameters
    ----------
    bam : str
          Path to bam file

    Returns
    -------
    counts : counter
             Counter with keys as number of times read maps
             and values as number of reads of that type
    """
    _create_bam_index(bam)
    bam = pysam.AlignmentFile(bam, 'rb')
    counts = Counter()
    for read in bam.fetch():
        if read.is_secondary:
            continue
        try:
            nh_count = Counter([dict(read.get_tags())['NH']])
        except KeyError:
            nh_count = Counter([1])
        counts += nh_count
    return counts


def metagene_coverage(bigwig,
                      htseq_f,
                      region_bed_f,
                      prefix=None,
                      master_offset=60,
                      top_n_meta=-1,
                      top_n_gene=10,
                      ignore_tx_version=True):
    bw = WigReader(bigwig)
    region_bed = pybedtools.BedTool(region_bed_f).sort().to_dataframe()

    # Group intervals by gene name
    cds_grouped = region_bed.groupby('name')

    # Get region sizes
    region_sizes = get_region_sizes(cds_grouped)
    ranked_genes = read_htseq(htseq_f, region_sizes, prefix)

    # Only consider genes which are in cds_grouped.keys
    ranked_genes = [gene for gene in ranked_genes if gene in cds_grouped.groups.keys()]
    if prefix:
        mkdir_p(os.path.dir(prefix))
        pickle.dump(ranked_genes,
                    open('{}_ranked_genes.pickle'.format(prefix), 'wb'),
                    pickle.HIGHEST_PROTOCOL)

    genewise_offsets = {}
    gene_position_counter = Counter()
    genewise_normalized_coverage = pd.Series()
    genewise_raw_coverage = pd.Series()

    if top_n_meta == -1:
        # Use all
        top_meta_genes = ranked_genes
    else:
        top_meta_genes = ranked_genes[:top_n_meta]
    topgene_normalized_coverage = pd.Series()
    topgene_position_counter = Counter()

    if top_n_gene == -1:
        # Use all genes! Not recommended
        top_genes = ranked_genes
    elif top_n_gene == 0:
        top_genes = []
    else:
        # Top  genes individual plot
        top_genes = ranked_genes[:top_n_gene]

    for gene_name, gene_group in cds_grouped:
        if ignore_tx_version:
            gene_name = re.sub(r'\.[0-9]+', '', gene_name)
        coverage_combined, gene_offset = gene_coverage(
            gene_group, bw, master_offset)

        # Generate individual plot for top genes
        if gene_name in top_genes:
            pickle.dump(coverage_combined,
                        open('{}_{}.pickle'.format(prefix, gene_name), 'wb'),
                        pickle.HIGHEST_PROTOCOL)

        # Generate top gene version metagene plot
        if gene_name in top_meta_genes:
            topgene_normalized_coverage = topgene_normalized_coverage.add(
                coverage_combined / coverage_combined.mean(), fill_value=0)

        genewise_normalized_coverage = genewise_normalized_coverage.add(
            coverage_combined / coverage_combined.mean(), fill_value=0)
        genewise_raw_coverage = genewise_raw_coverage.add(
            coverage_combined, fill_value=0)
        gene_position_counter += Counter(coverage_combined.index.tolist())
        genewise_offsets[gene_name] = gene_offset

    if len(gene_position_counter) != len(genewise_normalized_coverage):
        sys.exit(1)

    gene_position_counter = pd.Series(gene_position_counter)
    metagene_normalized_coverage = genewise_normalized_coverage.div(gene_position_counter)
    metagene_raw_coverage = genewise_raw_coverage
    pickle.dump(gene_position_counter,
                open('{}_gene_position_counter.pickle'.format(prefix), 'wb'),
                pickle.HIGHEST_PROTOCOL)
    pickle.dump(metagene_normalized_coverage,
                open('{}_metagene_normalized.pickle'.format(prefix), 'wb'),
                pickle.HIGHEST_PROTOCOL)
    pickle.dump(metagene_raw_coverage,
                open('{}_metagene_raw.pickle'.format(prefix), 'wb'),
                pickle.HIGHEST_PROTOCOL)
    pickle.dump(genewise_offsets,
                open('{}_genewise_offsets.pickle'.format(prefix), 'wb'),
                pickle.HIGHEST_PROTOCOL)

    if len(topgene_position_counter) != len(topgene_normalized_coverage):
        sys.exit(1)

    topgene_position_counter = pd.Series(topgene_position_counter)
    topgene_normalized_coverage = topgene_normalized_coverage.div(topgene_position_counter)

    pickle.dump(topgene_position_counter,
                open('{}_topgene_position_counter.pickle'.format(prefix), 'wb'),
                pickle.HIGHEST_PROTOCOL)
    pickle.dump(topgene_normalized_coverage,
                open('{}_topgene_normalized.pickle'.format(prefix), 'wb'),
                pickle.HIGHEST_PROTOCOL)
    return metagene_normalized_coverage

def read_htseq(htseq_f):
    """Read HTSeq file.

    Parameters
    ----------
    htseq_f : str
              Path to htseq counts file

    Returns
    -------
    htseq_df : dataframe
               HTseq counts as in a dataframe
    """
    htseq = pd.read_table(htseq_f, names=['name', 'counts']).set_index('name')
    htseq = htseq.iloc[:-5]
    if(htseq.shape[0] <= 10):
        print('Empty dataframe for : {}\n'.format(htseq_f))
        return None
    return htseq


def read_length_distribution(bam):
    """Count read lengths.

    Parameters
    ----------
    bam : str
          Path to bam file

    Returns
    -------
    lengths : counter
              Counter of read length and counts

    """
    _create_bam_index(bam)
    bam = pysam.AlignmentFile(bam, 'rb')
    return Counter([read.query_length for read in bam.fetch()
                    if _is_read_uniq_mapping(read)])


def summarize_counters(samplewise_dict):
    """Summarize gene counts for a collection of samples.

    Parameters
    ----------
    samplewise_dict : dict
                      A dictionary with key as sample name and value
                      as another dictionary of counts for each gene

    Returns
    -------
    totals : dict
             A dictionary with key as sample name and value as total gene count

    """
    totals = {}
    for key, sample_gene_dict in samplewise_dict.iteritems():
        totals[key] = np.nansum([np.nansum(d)
                                 for d in sample_gene_dict.values()])
    return totals


def unique_mapping_reads_count(bam):
    """Count number of mapped reads.

    Parameters
    ----------
    bam : str
          Path to bam file

    Returns
    -------
    n_mapped : int
               Count of mapped reads
    """
    _create_bam_index(bam)
    bam = pysam.AlignmentFile(bam, 'rb')
    n_mapped = 0
    for read in bam.fetch():
        if _is_read_uniq_mapping(read):
            n_mapped += 1
    return n_mapped
