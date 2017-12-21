"""Utilities for read counting operations."""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

from collections import Counter
from collections import OrderedDict
import subprocess
import sys

import numpy as np
import pandas as pd
import pybedtools
from pyfaidx import Fasta
import pysam
from scipy.stats import norm

from .wig import WigReader


def bedgraph_to_bigwig(bedgraph, chrom_sizes, bigwig):
    """Convert bedgraph to bigwig.

    Parameters
    ----------
    bedgraph : str
               Path to bedgraph file
    chrom_sizes : str
                  Path to genome chromosome sizes file
    bigwig : str
             Path to write bigwig file
    """
    cmds = ['bedSort', bedgraph, bedgraph]
    p = subprocess.Popen(cmds, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, universal_newlines=True)
    stdout, stderr = p.communicate()

    cmds = ['bedGraphToBigWig', bedgraph, chrom_sizes, bigwig]
    p = subprocess.Popen(cmds, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, universal_newlines=True)
    stdout, stderr = p.communicate()


def create_bedgraph(bam, strand='both', end_type='5prime', outfile=None):
    """Create bigwig from bam.

    Parameters
    ----------
    bam : str
          Path to bam file
    strand : str, optional
             Use fragments mapping to '+/-/both' strands
    end_type : str
               Use only end_type=5prime(5') or "3prime(3')"
    outfile : str, optional
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
    if outfile:
        with open(outfile, 'w') as outf:
            outf.write(str(genome_cov))
    return genome_cov


def fragment_enrichment(fragment_lengths, enrichment_range=range(28, 32)):
    """Calculate fragment enrichment for a certain range of lengths

    Parameters
    ----------
    fragment_lengths : Counter
                       A counter with fragment lengths and their counts
    enrichment_range : range
                       Range of fragments to concentrate upon

    Returns
    -------
    ratiro : float
             Enrichment in this range

    """
    if isinstance(fragment_lengths, Counter):
        fragment_lengths = pd.Series(fragment_lengths)
    rpf_signal = fragment_lengths[enrichment_range].sum()
    total_signal = fragment_lengths.sum()
    array = [[x] * y for x, y in sorted(fragment_lengths.iteritems())]
    mean_length, std_dev_length = norm.fit(
        np.concatenate(array).ravel().tolist())

    # mean_length_floor = np.floor(mean_length)
    # 1 - P(x1 < X <x2) = P(X<x1) + P(X>x2) = cdf(x1) + sf(x2)
    pvalue = norm.cdf(min(enrichment_range), mean_length, std_dev_length) + norm.sf(max(enrichment_range),
                                                                                    mean_length,
                                                                                    std_dev_length)
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
    intervals_for_fasta_query : list
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
    for index, coverage in enumerate(bw.query(intervals)):
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
    index_to_genomic_pos_map = pd.Series(coverage_combined.index.tolist(),
                                         index=np.arange(len(coverage_combined)) - gene_offset)
    intervals_for_fasta_query = []
    for pos in index_to_genomic_pos_map.values:
        intervals_for_fasta_query.append((chrom, pos, pos + 1, strand))
    coverage_combined = coverage_combined.reset_index(drop=True)
    coverage_combined = coverage_combined.rename(lambda x: x - gene_offset)

    return (coverage_combined, intervals_for_fasta_query,
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


def htseq_to_cpm(htseq_f, outfile=None):
    """Convert HTSeq counts to CPM.

    Parameters
    ----------
    htseq_f : str
              Path to HTseq counts file
    outfile : str, optional
              Path to output file
    Returns
    -------
    cpm : dataframe
          CPM

    """
    htseq = read_htseq(htseq_f)
    rate = htseq['counts']
    denom = rate.sum()
    cpm = rate/denom *1e6
    if outfile:
        pd.DataFrame(cpm, columns=['cpm']).to_csv(
            outfile, sep='\t', index=True, header=False)
    cpm = pd.DataFrame(cpm)#, columns=['cpm'])
    cpm.columns = ['cpm']
    cpm = cpm.sort_values(by='cpm', ascending=False)
    return cpm


def htseq_to_tpm(htseq_f, cds_bed_f, outfile=None):
    """Convert HTSeq counts to TPM.

    Parameters
    ----------
    htseq_f : str
              Path to HTseq counts file
    region_sizes : dict
                   Dict with keys as gene and values
                   as length (CDS/Exon) of that gene
    outfile : str, optional
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
    if outfile:
        pd.DataFrame(tpm, columns=['tpm']).to_csv(
            outfile, sep='\t', index=True, header=False)
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
    bam = pysam.AlignmentFile(bam, 'rb')
    counts = Counter()
    for query in bam.fetch():
        if query.is_secondary:
            continue
        try:
            nh_count = Counter([dict(query.get_tags())['NH']])
        except KeyError:
            nh_count = Counter([1])
        counts += nh_count
    return counts


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
    """Count fragment lengths.

    Parameters
    ----------
    bam : str
          Path to bam file

    Returns
    -------
    lengths : counter
              Counter of fragment length and counts

    """
    bam = pysam.AlignmentFile(bam, 'rb')
    return Counter([query.query_length for query in bam.fetch()
                    if query.mapping_quality == 255
                    and not query.is_secondary])


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
    bam = pysam.AlignmentFile(bam, 'rb')
    n_mapped = len([query for query in bam.fetch()
                    if query.mapping_quality == 255
                    and not query.is_secondary])
    return n_mapped
