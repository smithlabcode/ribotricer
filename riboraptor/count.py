"""Utilities for read counting operations."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from collections import Counter
from collections import OrderedDict
from contextlib import contextmanager
from functools import partial
import multiprocessing
import os
import pickle
import re
import subprocess
import sys

import HTSeq
import numpy as np
import pandas as pd
import pybedtools
from pyfaidx import Fasta
import pysam
import six
from scipy.stats import norm

from .genome import _get_sizes
from .genome import _get_bed
from .genome import __GENOMES_DB__

from .helpers import check_file_exists
from .helpers import collapse_bed_intervals
from .helpers import _fix_bed_coltype
from .helpers import list_to_ranges
from .helpers import load_pickle
from .helpers import mkdir_p
from .helpers import pad_or_truncate
from .helpers import pad_five_prime_or_truncate
from .helpers import summary_stats_two_arrays_welch

from .wig import WigReader

# Unmapped, Unmapped+Reverse strand, Not primary alignment,
# Not primary alignment + reverse strand, supplementary alignment

# Source: https://broadinstitute.github.io/picard/explain-flags.html
__SAM_NOT_UNIQ_FLAGS__ = [4, 20, 256, 272, 2048]
__PICKLE_PROTOCOL__ = 2
__LOOKUP_SUFFIX__ = '_genewise_lookup.pickle'


@contextmanager
def _poolcontext(*args, **kwargs):
    pool = multiprocessing.Pool(*args, **kwargs)
    yield pool
    pool.terminate()


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


def _bed_to_genomic_interval(bed):
    """ Converts bed file to genomic interval (htseq format) file

    Parameters
    ----------
    bed : str
          Path to bed file

    Returns
    -------
    interval : HTSeq.GenomicPosition
    """

    for interval in bed:
        yield HTSeq.GenomicPosition(
            str(interval.chrom), interval.start, str(interval.strand))


def get_closest(bam, regions, half_window_width=0):
    """For each read in bam find nearest region in bed.

    Parameters
    ----------
    bam : str
          Path to bam file
    regions : HTseq.GenomicPosition
              Genomic regions to get distance from
    half_window_width : int
                        Distance around region to record

    Returns
    -------
    profile : dataframe
    """
    profile = {}
    sorted_bam = HTSeq.BAM_Reader(bam)
    for x, tss in enumerate(_bed_to_genomic_interval(regions)):
        window = HTSeq.GenomicInterval(
            str(tss.chrom), tss.pos - half_window_width,
            tss.pos + half_window_width, str(tss.strand))
        for almnt in sorted_bam[window]:
            if almnt.iv.strand == '+':
                read_loc = almnt.iv.start - tss.pos
            else:
                read_loc = tss.pos - almnt.iv.end
            length = sum(cigar.size for cigar in almnt.cigar
                         if cigar.type == 'M')
            profile[almnt.read.name] = {'dist': read_loc, 'length': length}
    return pd.DataFrame(profile)


def gene_coverage(gene_name,
                  bed,
                  bw,
                  gene_group=None,
                  offset_5p=0,
                  offset_3p=0,
                  collapse_intervals=True):
    """Get gene coverage.

    Parameters
    ----------
    gene_name : str
                Gene name
    bed : str
          Path to CDS or 5'UTR or 3'UTR bed
    bw : str
         Path to bigwig to fetch the scores from
    offset_5p : int (positive)
                Number of bases to count upstream (5')
    offset_3p : int (positive)
                Number of bases to count downstream (3')
    collapse_intervals : bool
                         Should bed be collapsed based on gene name

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
    assert (offset_5p >= 0)
    assert (offset_3p >= 0)
    if not isinstance(bw, WigReader):
        bw = WigReader(bw)
    chromosome_lengths = bw.get_chromosomes
    bed_df = pybedtools.BedTool(bed).to_dataframe()
    bed_df = _fix_bed_coltype(bed_df)
    assert gene_name in bed_df['name'].tolist()
    if gene_group is None:
        gene_group = bed_df[bed_df['name'] == gene_name]

    if len(gene_group['strand'].unique()) != 1:
        raise RuntimeError('Multiple strands?: {}'.format(gene_group))
    if len(gene_group['chrom'].unique()) != 1:
        raise RuntimeError('Chrom not unique for: {}'.format(gene_group))
    chrom = gene_group['chrom'].unique()[0]
    strand = gene_group['strand'].unique()[0]

    # Collect all intervals at once
    intervals = list(
        zip(gene_group['chrom'], gene_group['start'], gene_group['end'],
            gene_group['strand']))

    query_intervals, fasta_onebased_intervals, gene_offset_5p, gene_offset_3p = collapse_bed_intervals(
        intervals, chromosome_lengths, offset_5p, offset_3p)
    try:
        interval_lookup_pickle = bed + __LOOKUP_SUFFIX__
        if check_file_exists(interval_lookup_pickle):
            query_intervals = load_pickle(interval_lookup_pickle)[gene_name]
            fasta_onebased_intervals = []
            for interval in query_intervals:
                # 1-based ending for both
                # Original bwas based on range and hence was 0-based
                fasta_onebased_intervals.append((interval[0], interval[1] + 1,
                                                 interval[2] + 1, interval[3]))
        elif collapse_intervals:
            query_intervals, fasta_onebased_intervals, gene_offset_5p, gene_offset_3p = collapse_bed_intervals(
                intervals, chromosome_lengths, offset_5p, offset_3p)
        else:
            query_intervals = intervals
    except Exception as error:
        sys.stderr.write(
            'Error reading {} \t {} : {}'.format(gene_name, intervals, error))
        sys.exit(1)
    coverage = bw.query(query_intervals)
    if len(coverage) == 0:
        return (pd.Series([]), pd.Series([]), pd.Series([]), 0, 0)
    coverage_combined = list(coverage[0])
    for coverage in coverage[1:]:
        coverage_combined += list(coverage)
    coverage_combined = np.array(coverage_combined).flatten()
    coverage_combined = pd.Series(
        coverage_combined,
        index=np.arange(-gene_offset_5p,
                        len(coverage_combined) - gene_offset_5p))
    return coverage_combined, fasta_onebased_intervals, gene_offset_5p, gene_offset_3p


def gene_coverage_sum(gene_name, bed, bw, collapse_intervals=True):
    """Keep track of only the sum

    Parameters
    ----------
    gene_name : str
                Name of gene
    bed : str
          Path to bed file
    bw : str
         Path to bigwig file
    collapse_intervals : bool
                         Should the intervals be collapsed based on the 'name' column in gene
                         This should be set to False for things like tRNA
                         where the tRNA can span multiple chromosomes

    """
    gene_group = None
    if not isinstance(bw, WigReader):
        bw = WigReader(bw)
    chromsome_lengths = bw.get_chromosomes
    if not isinstance(bed, pd.DataFrame):
        bed = pybedtools.BedTool(bed).to_dataframe()
        bed = _fix_bed_coltype(bed)
    assert gene_name in bed['name'].tolist()

    if gene_group is None:
        gene_group = bed[bed['name'] == gene_name]

    chrom = gene_group['chrom'].unique()[0]
    strand = gene_group['strand'].unique()[0]

    intervals = list(
        zip(gene_group['chrom'], gene_group['start'], gene_group['end'],
            gene_group['strand']))

    try:
        interval_lookup_pickle = bed + __LOOKUP_SUFFIX__
        if check_file_exists(interval_lookup_pickle):
            query_intervals = load_pickle(interval_lookup_pickle)[gene_name]
        elif collapse_intervals:
            query_intervals, fasta_onebased_intervals, gene_offset_5p, gene_offset_3p = collapse_bed_intervals(
                intervals)
        else:
            query_intervals = intervals
    except Exception as error:
        sys.stderr.write(
            'Error reading {} \t {} : {}'.format(gene_name, intervals, error))
        sys.exit(1)
    coverage = bw.query(query_intervals)
    if len(coverage) == 0:
        ## This gene does not exist in the bigwig
        return (np.nan, sum([x[2] - x[1] for x in query_intervals]))
    coverage_combined = list(coverage[0])
    for coverage in coverage[1:]:
        coverage_combined += list(coverage)
    coverage_combined = np.array(coverage_combined).flatten()
    coverage_combined = pd.Series(
        coverage_combined,
        index=np.arange(-gene_offset_5p,
                        len(coverage_combined) - gene_offset_5p))
    #print (coverage_combined.sum(skipna=true))
    return coverage_combined.sum(skipna=True), len(coverage_combined)


def pickle_bed_file(bed, collapse_intervals=True):
    """Create a lookup pickle file for genewise CDS/UTR coordinates.

    In order to prevent recalculating the coordinates that should be fetched
    for each genes' CDS or UTR regions, they can be stored in a pickle file.

    Parameters
    ----------
    bed : string
          Path to bed file
    collapse_intervals : bool
                         Should the intervals be collapsed based on the 'name' column in gene
                         This should be set to False for things like tRNA
                         where the tRNA can span multiple chromosomes
    """
    bed_df = pybedtools.BedTool(bed).sort().to_dataframe()
    bed_df = _fix_bed_coltype(bed_df)
    bed_grouped = bed_df.groupby('name')
    genewise_collapsed_intervals = {}
    for gene_name, bed_group in bed_grouped:
        intervals = list(
            zip(bed_group['chrom'], bed_group['start'], bed_group['end'],
                bed_group['strand']))
        if collapse_intervals:
            intervals, _, _, _ = collapse_bed_intervals(intervals)
        genewise_collapsed_intervals[gene_name] = intervals
    outfile = bed + __LOOKUP_SUFFIX__
    pickle.dump(genewise_collapsed_intervals, open(outfile, 'wb'))
    return genewise_collapsed_intervals


def count_reads_per_gene(bw,
                         bed,
                         prefix=None,
                         n_cores=16,
                         collapse_intervals=True):
    """Count number of reads following in each region.

    Parameters
    ----------
    bw : str
          Path to bigWig file
    bed : pybedtools.BedTool or str
          Genomic regions to get distance from
    prefix : str
             Prefix to output pickle files
    n_cores : int
              Use multiple cores (Default: 16). Set to 1 to disable multiprocessing
    collapse_intervals : bool
                         Should the intervals be collapsed based on the 'name' column in gene
                         This should be set to False for things like tRNA
                         where the tRNA can span multiple chromosomes

    Returns
    -------

    counts_by_region : Series
                       Series with counts indexed by gene id
    region_lengths : Series
                     Series with gene lengths
    counts_normalized_by_length : Series
                                  Series with normalized counts
    """
    counts_by_region = OrderedDict()
    length_by_region = OrderedDict()
    if isinstance(bed, six.string_types):
        bed = pybedtools.BedTool(bed).sort()
    bed_df = bed.to_dataframe()
    bed_df = _fix_bed_coltype(bed_df)
    gene_names = list(sorted(bed_df['name'].unique()))
    if n_cores == 1:
        for gene_name in gene_names:
            coverage_combined, _, _, _, _ = gene_coverage(
                gene_name, bed, bw, collapse_intervals)
            length = len(coverage_combined)
            counts_by_region[gene_name] = coverage_combined.sum()
            length_by_region[gene_name] = length
    else:
        with _poolcontext(processes=n_cores) as pool:
            results = pool.map(
                partial(
                    gene_coverage_sum,
                    bed=bed,
                    bw=bw,
                    collapse_intervals=collapse_intervals), list(gene_names))
        counts_by_region = OrderedDict(
            list(zip(list(gene_names), list([x[0] for x in results]))))
        length_by_region = OrderedDict(
            list(zip(list(gene_names), list([x[1] for x in results]))))
    counts_by_region = pd.Series(counts_by_region)
    length_by_region = pd.Series(length_by_region)
    counts_normalized_by_length = counts_by_region.div(length_by_region)
    if prefix:
        mkdir_p(os.path.dirname(prefix))
        df = pd.concat(
            [counts_by_region, length_by_region, counts_normalized_by_length],
            axis=1)
        df.columns = ['counts', 'length', 'normalized_counts']
        df.to_csv(
            '{}.tsv'.format(prefix), index=True, header=True, sep=str('\t'))
        pickle.dump(counts_by_region,
                    open('{}_counts.pickle'.format(prefix), 'wb'),
                    __PICKLE_PROTOCOL__)
        pickle.dump(length_by_region,
                    open('{}_lengths.pickle'.format(prefix), 'wb'),
                    __PICKLE_PROTOCOL__)
        pickle.dump(counts_normalized_by_length,
                    open('{}_counts_lengths_normalized.pickle'.format(prefix),
                         'wb'), __PICKLE_PROTOCOL__)
    return counts_by_region, length_by_region, counts_normalized_by_length


def bedgraph_to_bigwig(bedgraph, sizes, saveto, input_is_stream=False):
    """Convert bedgraph to bigwig.

    Parameters
    ----------
    bedgraph : str
               Path to bedgraph file
    sizes : str
            Path to genome chromosome sizes file
            or genome name
    saveto : str
             Path to write bigwig file
    input_is_stream : bool
                      True if input is sent through stdin
    """
    if input_is_stream:
        total_lines = len(bedgraph)
        with open(os.path.splitext(saveto)[0] + '.bg', 'w') as fp:
            for index, line in enumerate(bedgraph):
                if index == (total_lines - 1):
                    fp.write(line.rstrip())
                else:
                    fp.write(line)
            filename = str(fp.name)
        bedgraph = filename

    if not os.path.isfile(sizes):
        if sizes in __GENOMES_DB__:
            sizes = _get_sizes(sizes)
        else:
            raise RuntimeError('Could not load size for {}'.format(sizes))

    cmds = ['bedSort', bedgraph, bedgraph]
    p = subprocess.Popen(
        cmds,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True)
    stdout, stderr = p.communicate()
    rc = p.returncode
    if rc != 0:
        raise RuntimeError(
            'Error running bedSort.\nstdout : {} \n stderr : {}'.format(
                stdout, stderr))

    cmds = ['bedGraphToBigWig', bedgraph, sizes, saveto]
    p = subprocess.Popen(
        cmds,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True)
    stdout, stderr = p.communicate()
    rc = p.returncode
    if rc != 0:
        raise RuntimeError(
            'Error running bedSort.\nstdout : {} \n stderr : {}'.format(
                stdout, stderr))


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
        genome_cov = bed.genome_coverage(
            bg=True, strand=strand, additional_args=extra_args)
    else:
        genome_cov = bed.genome_coverage(bg=True, additional_args=extra_args)
    if saveto:
        with open(saveto, 'w') as outf:
            outf.write(str(genome_cov))
    return str(genome_cov)


def count_reads_in_features(feature_bed,
                            bam,
                            force_strandedness=False,
                            use_multiprocessing=False):
    """Count reads overlapping features.

    Parameters
    ----------
    feature_bed : str
                   Path to features bed file
    bam : str
          Path to bam file
    force_strandedness : bool
                         Should count feature only if on the same strand
    use_multiprocessing : bool
                          True if multiprocessing mode
    Returns
    -------
    counts : int
             Number of intersection between bam and bed
    """
    if not isinstance(bam, pybedtools.BedTool):
        bam = pybedtools.BedTool(bam)
    if use_multiprocessing:
        if force_strandedness:
            count = bam.intersect(
                b=feature_bed, bed=True, stream=True,
                additional_args='-s').count()
        else:
            count = bam.intersect(b=feature_bed, bed=True, stream=True).count()

    else:
        if not isinstance(feature_bed, pybedtools.BedTool):
            feature_bed = pybedtools.BedTool(feature_bed)
        if force_strandedness:
            count = bam.intersect(b=feature_bed, additional_args='-s').count()
        else:
            count = bam.intersect(b=feature_bed, stream=False).count()
    return count


def count_feature_genewise(feature_bed,
                           bam,
                           force_strandedness=False,
                           use_multiprocessing=False):
    """Count features genewise.

    Parameters
    ----------
    bam : str
          Path to bam file
    feature_bed : str
                   Path to features bed file

    Returns
    -------
    counts : dict
             Genewise feature counts
    """
    if not isinstance(feature_bed, pybedtools.BedTool):
        feature_bed = pybedtools.BedTool()
    feature_bed_df = feature_bed.to_dataframe()
    feature_bed_df = _fix_bed_coltype(feature_bed_df)

    feature_bed_df_grouped = feature_bed_df.groupby('name')
    genewise_beds = OrderedDict()
    for gene_name, gene_group in feature_bed_df_grouped:
        gene_bed = pybedtools.BedTool.from_dataframe(
            gene_group).sort().saveas().fn
        genewise_beds[gene_name] = gene_bed
    """
    with _poolcontext(processes=3) as pool:
        results = pool.map(partial(count_reads_in_features,
                                   bam=bam,
                                   force_strandedness=force_strandedness,
                                   use_multiprocessing=True), list(genewise_beds.values()))
    counts = OrderedDict(zip(list(gene_wise_beds.keys()), results))
    """
    counts = OrderedDict()
    for gene_name, gene_bed in six.iteritems(genewise_beds):
        counts[gene_name] = count_reads_in_features(
            pybedtools.BedTool(gene_bed), bam, force_strandedness, False)
    return counts


def count_utr5_utr3_cds(bam,
                        utr5_bed=None,
                        cds_bed=None,
                        utr3_bed=None,
                        genome=None,
                        force_strandedness=False,
                        genewise=False,
                        saveto=None,
                        use_multiprocessing=False):
    """One shot counts over UTR5/UTR3/CDS.

    Parameters
    ----------
    bam : str
          Path to bam file
    utr5_bed : str
               Path to 5'UTR feature bed file
    utr3_bed : str
               Path to 3'UTR feature bed file
    cds_bed : str
              Path to CDS feature bed file
    saveto : str, optional
             Path to output file
    use_multiprocessing : bool
                          SHould use multiprocessing?
                          Not been well tested if it really helps

    Returns
    -------
    counts : dict
             Dict with keys as feature type and counts as values

    """
    bam = pybedtools.BedTool(bam)
    order = ['utr5', 'cds', 'utr3']
    if genome:
        utr5_bed = _get_bed('utr5', genome)
        cds_bed = _get_bed('cds', genome)
        utr3_bed = _get_bed('utr3', genome)
    feature_beds = [utr5_bed, cds_bed, utr3_bed]
    for index, bed in enumerate(feature_beds):
        if bed is None:
            raise RuntimeError('{} cannot be none.'.format(order[index]))
    if genewise:
        counts = OrderedDict()
        for index, bed in enumerate(feature_beds):
            counts[order[index]] = count_feature_genewise(
                pybedtools.BedTool(bed), bam, force_strandedness)
    else:
        if use_multiprocessing:
            with _poolcontext(processes=3) as pool:
                results = pool.map(
                    partial(
                        count_reads_in_features,
                        bam=bam,
                        force_strandedness=force_strandedness,
                        use_multiprocessing=True), feature_beds)
            counts = OrderedDict(list(zip(order, results)))
        else:
            counts = OrderedDict()
            for index, bed in enumerate(feature_beds):
                counts[order[index]] = count_reads_in_features(
                    pybedtools.BedTool(bed), bam, force_strandedness)
    if saveto:
        mkdir_p(os.path.dirname(saveto))
        pickle.dump(counts,
                    open('{}.pickle'.format(saveto), 'wb'),
                    __PICKLE_PROTOCOL__)
    return counts


def diff_region_enrichment(numerator, denominator, prefix):
    """Calculate enrichment of counts of one region over another.

    Parameters
    ----------
    numerator : str
                Path to pickle file
    denominator : str
                  Path to pickle file
    prefix : str
             Prefix to save pickles to

    Returns
    --------
    enrichment : series
    """
    numerator = load_pickle(numerator)
    denominator = load_pickle(denominator)

    enrichment = numerator.divide(denominator)
    if prefix:
        mkdir_p(os.path.dirname(prefix))
        pickle.dump(enrichment,
                    open('{}.pickle'.format(prefix), 'wb'),
                    __PICKLE_PROTOCOL__)
    return enrichment


def read_enrichment(read_lengths,
                    enrichment_range=list(range(28, 33)),
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
            Enrichment in this range (Scae 0-1)

    """
    if input_is_file:
        if not check_file_exists(read_lengths):
            raise RuntimeError('{} does not exist.'.format(read_lengths))
        try:
            read_lengths = load_pickle(read_lengths)
        except KeyError:
            read_lengths = pd.read_table(
                read_lengths, names=['frag_len', 'frag_count'], sep='\t')
            read_lengths = pd.Series(
                read_lengths.frag_count.tolist(),
                index=read_lengths.frag_len.tolist())
    elif input_is_stream:
        counter = {}
        for line in read_lengths:
            splitted = list([int(x) for x in line.strip().split('\t')])
            counter[splitted[0]] = splitted[1]
        read_lengths = Counter(counter)
    if isinstance(read_lengths, Counter):
        read_lengths = pd.Series(read_lengths)
    if isinstance(enrichment_range, six.string_types):
        splitted = list([int(x) for x in enrichment_range.strip().split('-')])
        enrichment_range = list(range(splitted[0], splitted[1] + 1))
    rpf_signal = read_lengths[enrichment_range].sum()
    total_signal = read_lengths.sum()
    array = [[x] * y for x, y in list(sorted(six.iteritems(read_lengths)))]
    mean_length, std_dev_length = norm.fit(
        np.concatenate(array).ravel().tolist())

    # mean_length_floor = np.floor(mean_length)
    # 1 - P(x1 < X <x2) = P(X<x1) + P(X>x2) = cdf(x1) + sf(x2)
    cdf_min = norm.cdf(min(enrichment_range), mean_length, std_dev_length)
    sf_max = norm.sf(max(enrichment_range), mean_length, std_dev_length)
    pvalue = cdf_min + sf_max
    ratio = rpf_signal / float(total_signal)
    return ratio, pvalue


def interval_coverage(bw, intervals):
    """Get coverage at custom intervals

    Parameters
    ----------
    bw : str
         Path to bigwig file
    intervals : list of tuples
                [(chrom, start, stop, strand)]

    Returns
    -------
    coverage : list of series
               Coverage for each interval, so that
               it is sorted oritentation wise
    """
    if not isinstance(bw, WigReader):
        bw = WigReader(bw)
    if isinstance(intervals, tuple):
        # Only empty interval
        intervals = [intervals]
    interval_coverage_list = []
    chrom = intervals[0][0]
    for index, coverage in enumerate(bw.query(intervals)):
        strand = intervals[index][3]
        if strand == '+':
            series_range = list(
                range(intervals[index][1], intervals[index][2]))
        elif strand == '-':
            series_range = list(
                range(intervals[index][2], intervals[index][1], -1))
        series = pd.Series(coverage, index=series_range)
        interval_coverage_list.append(series)
    coverage_combined = interval_coverage_list[0]
    for interval_coverage in interval_coverage_list[1:]:
        coverage_combined = coverage_combined.combine_first(interval_coverage)
    coverage_combined = coverage_combined.fillna(0)
    coverage_index = np.arange(len(coverage_combined))
    index_to_genomic_pos_map = pd.Series(
        coverage_combined.index.tolist(), index=coverage_index)
    intervals_for_fasta_read = []
    for pos in index_to_genomic_pos_map.values:
        intervals_for_fasta_read.append((chrom, pos, pos + 1, strand))
    if strand == '+':
        coverage_combined = coverage_combined.reset_index(drop=True)
    else:
        coverage_combined = coverage_combined.sort_index(
            ascending=False).reset_index(drop=True)

    return (coverage_combined, intervals_for_fasta_read,
            index_to_genomic_pos_map)


def export_gene_coverages(bigwig,
                          region_bed_f,
                          prefix,
                          offset_5p=60,
                          offset_3p=0,
                          ignore_tx_version=True):
    """Export all gene coverages.

    Parameters
    ----------
    bigwig : str
             Path to bigwig file
    region_bed_f : str
                   Path to region bed file (CDS/3'UTR/5'UTR)
                   with bed name column as gene
    prefix : str
             Prefix to write output file
    offset_5p : int
             number of bases to count upstream (5')
    offset_30 : int
                number of bases to count downstream (3')
    ignore_tx_version : bool
                        Should versions be ignored for gene names

    Returns
    -------
    gene_profiles: file
                   with the following format:
                   gene1\t5poffset1\t3poffset1\tlength1\tmean1\tmedian1\tstdev1\tcnt1_1 cnt1_2 cnt1_3 ...\n
                   gene2\t5poffset2\t3poffset2\tlength2\tmean2\tmedian2\tstdev2\tcnt2_1 cnt2_2 cnt2_3 cnt2_4 ...\n
    """
    bw = WigReader(bigwig)
    if region_bed_f.lower().split('_')[0] in __GENOMES_DB__:

        genome, region_type = region_bed_f.lower().split('_')
        region_bed_f = _get_bed(region_type, genome)

    region_bed = pybedtools.BedTool(region_bed_f).sort().to_dataframe()
    region_bed['chrom'] = region_bed['chrom'].astype(str)
    region_bed['name'] = region_bed['name'].astype(str)
    # Group intervals by gene name
    cds_grouped = region_bed.groupby('name')

    with open('{}_gene_coverages.tsv'.format(prefix), 'w') as outfile:
        outfile.write(
            'gene_name\toffset_5p\toffset_3p\tlength\tmean\tmedian\tstdev\tcount\n'
        )
        for gene_name, gene_group in cds_grouped:
            if ignore_tx_version:
                gene_name = re.sub(r'\.[0-9]+', '', gene_name)
            gene_cov, _, gene_offset_5p, gene_offset_3p = gene_coverage(
                gene_name, region_bed_f, bw, gene_group, offset_5p, offset_3p)
            gene_cov = gene_cov.fillna(0)
            count = gene_cov.tolist()
            length = len(count)
            count_array = np.array(count)
            mean = np.nanmean(count_array)
            median = np.nanmedian(count_array)
            stdev = np.nanstd(count_array)
            # Write olny useful etnries which has mean > 0
            if mean > 0:
                outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    gene_name, gene_offset_5p,
                    int(gene_offset_3p), length, mean, median, stdev, count))


def export_single_gene_coverage(bigwig,
                                region_bed_f,
                                gene_name,
                                prefix,
                                offset_5p=60,
                                offset_3p=0,
                                ignore_tx_version=True):
    """Export single gene coverage to tsv

    Parameters
    ----------
    bigwig : str
             Path to bigwig file
    region_bed_f : str
                   Path to region bed file (CDS/3'UTR/5'UTR)
                   with bed name column as gene
    prefix : str
             Prefix to write output file
    offset : int
             Number of bases to count upstream
    ignore_tx_version : bool
                        Should versions be ignored for gene names

    Returns
    -------
    gene_profile: tsv with first column as index, next coumn as value
    """
    bw = WigReader(bigwig)
    if region_bed_f.lower().split('_')[0] in __GENOMES_DB__:

        genome, region_type = region_bed_f.lower().split('_')
        region_bed_f = _get_bed(region_type, genome)

    region_bed = pybedtools.BedTool(region_bed_f).sort().to_dataframe()
    region_bed['chrom'] = region_bed['chrom'].astype(str)
    region_bed['name'] = region_bed['name'].astype(str)
    # Group intervals by gene name
    cds_grouped = region_bed.groupby('name')
    if ignore_tx_version:
        gene_name = re.sub(r'\.[0-9]+', '', gene_name)
    if gene_name not in region_bed['name'].tolist():
        print(region_bed['name'])
        raise RuntimeError('{} not in bed'.format(gene_name))

    gene_group = cds_grouped.get_group(gene_name)
    gene_cov, _, _, gene_offset_5p, gene_offset_3p = gene_coverage(
        gene_name, region_bed, bw, gene_group, offset_5p, offset_3p)
    gene_cov.to_csv(
        path='{}_coverage.tsv'.format(prefix), sep=str('\t'), index=True)


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
    bed['chrom'] = bed['chrom'].astype(str)
    bed['name'] = bed['name'].astype(str)
    region_bed_grouped = bed.groupby('name')
    region_sizes = {}
    for gene_name, gene_group in region_bed_grouped:
        # Get rid of trailing dots
        # TODO: this can be letha for genomes with names
        # such as orf19.123
        # gene_name = re.sub(r'\.[0-9]+', '', gene_name)
        # Collect all intervals at once
        intervals = list(
            zip(gene_group['chrom'], gene_group['start'], gene_group['end'],
                gene_group['strand']))
        collapsed_intervals, _, _, _ = collapse_bed_intervals(intervals)
        for interval in collapsed_intervals:
            if gene_name not in region_sizes:
                # End is always 1-based so does not require +1
                region_sizes[gene_name] = interval[2] - interval[1]
            else:
                region_sizes[gene_name] += interval[2] - interval[1]
    sizes_df = pd.DataFrame.from_dict(
        region_sizes, orient='index').reset_index()
    sizes_df.columns = ['name', 'length']
    sizes_df.sort_values(by='length', ascending=False, inplace=True)
    region_sizes = OrderedDict(
        [tuple(x) for x in sizes_df[['name', 'length']].values])
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
        pd.DataFrame(
            cpm, columns=['cpm']).to_csv(
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
    if cds_bed_f.lower().split('_')[0] in __GENOMES_DB__:
        cds_bed_f = _get_bed('cds', cds_bed_f)
    cds_bed_sizes = get_region_sizes(cds_bed_f)
    htseq = read_htseq(htseq_f)
    rate = np.log(htseq['counts']).subtract(np.log(pd.Series(cds_bed_sizes)))
    denom = np.log(np.sum(np.exp(rate)))
    tpm = np.exp(rate - denom + np.log(1e6))
    tpm = pd.DataFrame(tpm, columns=['tpm'])
    tpm = tpm.sort_values(by='tpm', ascending=False)
    if saveto:
        tpm.to_csv(saveto, sep=str('\t'), index=True, header=False)
    return tpm


def mapping_reads_summary(bam, prefix):
    """Count number of mapped reads.

    Parameters
    ----------
    bam : str
          Path to bam file
    prefix : str
             Prefix to save pickle to (optional)

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
    if prefix:
        mkdir_p(os.path.dirname(prefix))
        pickle.dump(counts,
                    open('{}.pickle'.format(prefix), 'wb'),
                    __PICKLE_PROTOCOL__)
    return counts


def metagene_coverage(bigwig,
                      region_bed_f,
                      max_positions=None,
                      htseq_f=None,
                      prefix=None,
                      offset_5p=60,
                      offset_3p=0,
                      top_n_meta=-1,
                      top_n_gene=10,
                      ignore_tx_version=True,
                      calculate_welch=False):
    """Calculate metagene coverage.

    Parameters
    ----------
    bigwig : str
             Path to bigwig file
    region_bed_f : str
                   Path to region bed file (CDS/3'UTR/5'UTR)
                   with bed name column as gene
                   or a genome name (hg38_cds, hg38_utr3, hg38_utr5)
    max_positions: int
                   Number of positions to consider while
                   calculating the normalized coverage
                   -1 for the entore length of gene
                   Higher values lead to slower implementation
    htseq_f : str
              Path to htseq-counts file
    prefix : str
             Prefix to write output files
    offset_5p : int
                Number of bases to offset upstream(5')
    offset_3p : int
                Number of bases to offset downstream(3')
    top_n_meta : int
                 Total number of top expressed genes
                 to use for calculating metagene profile
    top_n_gene : int
                 Number of gene profiles to output
    ignore_tx_version : bool
                 Should versions be ignored for gene names

    Returns
    -------
    metagene_profile : series
                       Metagene profile
    """
    welch_var = None
    bw = WigReader(bigwig)
    if region_bed_f.lower().split('_')[0] in __GENOMES_DB__:

        genome, region_type = region_bed_f.lower().split('_')
        region_bed_f = _get_bed(region_type, genome)

    region_bed = pybedtools.BedTool(region_bed_f).sort().to_dataframe()
    region_bed['chrom'] = region_bed['chrom'].astype(str)
    region_bed['name'] = region_bed['name'].astype(str)
    # Group intervals by gene name
    cds_grouped = region_bed.groupby('name')

    # Get region sizes
    if htseq_f:
        ranked_genes = htseq_to_tpm(htseq_f, region_bed_f).index.tolist()
        # region_sizes = get_region_sizes(region_bed_f)
        # Only consider genes which are in cds_grouped.keys
        ranked_genes = [
            gene for gene in ranked_genes
            if gene in list(cds_grouped.groups.keys())
        ]
    else:
        # Use all genes when no htseq present
        ranked_genes = list(cds_grouped.groups.keys())
    genewise_offsets = {}
    gene_position_counter = Counter()
    genewise_normalized_coverage = pd.Series()
    genewise_raw_coverage = pd.Series()

    if top_n_meta == -1:
        # Use all
        top_meta_genes = ranked_genes
    else:
        # Not useful if htseq-counts missing
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
    index = 0
    for gene_name, gene_group in cds_grouped:
        if ignore_tx_version:
            gene_name = re.sub(r'\.[0-9]+', '', gene_name)
        gene_cov, _, gene_offset_5p, gene_offset_3p = gene_coverage(
            gene_name, region_bed_f, bw, gene_group, offset_5p, offset_3p)
        if max_positions is not None and len(gene_cov.index) > 0:
            min_index = min(gene_cov.index.tolist())
            gene_length = max(gene_cov.index.tolist())
            # gene_length = len(gene_cov.inex) + min_index

            # Keep only important positions
            gene_cov = gene_cov[np.arange(min_index,
                                          min(gene_length, max_positions))]

        # Generate individual plot for top genes
        if gene_name in top_genes and prefix:
            mkdir_p(os.path.dirname(prefix))
            pickle.dump(gene_cov,
                        open('{}_{}.pickle'.format(prefix, gene_name), 'wb'),
                        __PICKLE_PROTOCOL__)

        # Generate top gene version metagene plot
        norm_cov = gene_cov / gene_cov.mean()
        if index == 0:
            previous_gene_coverage = norm_cov
            previous_variance = None
            old_n_counter = Counter(gene_cov.index.tolist())
            welch_carried_forward = None
        if index == 1:
            old_n_counter += Counter(gene_cov.index.tolist())

        if index >= 1:
            new_gene_coverage = norm_cov
        if gene_name in top_meta_genes:
            topgene_normalized_coverage = topgene_normalized_coverage.add(
                norm_cov, fill_value=0)
            topgene_position_counter += Counter(gene_cov.index.tolist())

        genewise_normalized_coverage = genewise_normalized_coverage.add(
            norm_cov, fill_value=0)
        if index > 1 and len(new_gene_coverage) > 1 and len(
                previous_gene_coverage) > 1 and calculate_welch:
            welch_mean, welch_var, welch_n_obs, welch_carried_forward = summary_stats_two_arrays_welch(
                old_mean_array=previous_gene_coverage,
                new_array=new_gene_coverage,
                old_var_array=previous_variance,
                old_n_counter=old_n_counter,
                carried_forward_observations=welch_carried_forward)
            previous_gene_coverage = welch_mean.copy()
            previous_variance = welch_var.copy()
            old_n_counter = welch_n_obs.copy()

        genewise_raw_coverage = genewise_raw_coverage.add(
            gene_cov, fill_value=0)
        gene_position_counter += Counter(gene_cov.index.tolist())
        """
        if index > 1:
            if not (pd.Series(welch_n_obs) == pd.Series(gene_position_counter)
                    ).all():
                print('welch_n_ob n obs :{}'.format(pd.Series(welch_n_obs)))
                print('gene pos obs :{}'.format(
                    pd.Series(gene_position_counter)))
        """
        genewise_offsets[gene_name] = gene_offset_5p
        index += 1

    if len(gene_position_counter) != len(genewise_normalized_coverage):
        raise RuntimeError('Gene normalizaed counter mismatch')
        sys.exit(1)

    gene_position_counter = pd.Series(gene_position_counter)
    metagene_normalized_coverage = genewise_normalized_coverage.div(
        gene_position_counter)
    metagene_raw_coverage = genewise_raw_coverage

    if len(topgene_position_counter) != len(topgene_normalized_coverage):
        raise RuntimeError('Top gene normalizaed counter mismatch')
        sys.exit(1)

    topgene_position_counter = pd.Series(topgene_position_counter)
    topgene_normalized_coverage = topgene_normalized_coverage.div(
        topgene_position_counter)

    if prefix:
        mkdir_p(os.path.dirname(prefix))
        pickle.dump(ranked_genes,
                    open('{}_ranked_genes.pickle'.format(prefix), 'wb'),
                    __PICKLE_PROTOCOL__)
        pickle.dump(gene_position_counter,
                    open('{}_gene_position_counter.pickle'.format(prefix),
                         'wb'), __PICKLE_PROTOCOL__)
        pickle.dump(metagene_normalized_coverage,
                    open('{}_metagene_normalized.pickle'.format(prefix), 'wb'),
                    __PICKLE_PROTOCOL__)
        pickle.dump(metagene_raw_coverage,
                    open('{}_metagene_raw.pickle'.format(prefix), 'wb'),
                    __PICKLE_PROTOCOL__)
        if welch_var is not None:
            pickle.dump(welch_var,
                        open('{}_metagene_var.pickle'.format(prefix), 'wb'),
                        __PICKLE_PROTOCOL__)
        pickle.dump(genewise_offsets,
                    open('{}_genewise_offsets.pickle'.format(prefix), 'wb'),
                    __PICKLE_PROTOCOL__)
        pickle.dump(topgene_position_counter,
                    open('{}_topgene_position_counter.pickle'.format(prefix),
                         'wb'), __PICKLE_PROTOCOL__)
        pickle.dump(topgene_normalized_coverage,
                    open('{}_topgene_normalized.pickle'.format(prefix), 'wb'),
                    __PICKLE_PROTOCOL__)
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
    htseq = pd.read_table(htseq_f, names=['name', 'counts'])
    htseq['name'] = htseq['name'].astype(str)
    htseq = htseq.set_index('name')
    htseq = htseq.iloc[:-5]
    if (htseq.shape[0] <= 10):
        sys.stderr.write('Empty dataframe for : {}\n'.format(htseq_f))
        return None
    return htseq


def read_length_distribution(bam, prefix):
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
    read_counts = Counter([
        read.query_length for read in bam.fetch()
        if _is_read_uniq_mapping(read)
    ])
    if prefix:
        mkdir_p(os.path.dirname(prefix))
        pickle.dump(
            pd.Series(read_counts),
            open('{}.pickle'.format(prefix), 'wb'), __PICKLE_PROTOCOL__)

    return read_counts


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
    for key, sample_dict in six.iteritems(samplewise_dict):
        totals[key] = np.nansum(
            [np.nansum(d) for d in list(sample_dict.values)])
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


def collapse_gene_coverage_to_metagene(gene_coverages,
                                       target_length,
                                       outfile=None):
    """Collapse gene coverages to specific target length.

    Parameters
    ----------
    gene_coverages : string
                     Path to gene coverages.tsv
    target_lenght : int
                    Collapse to target length

    Returns
    -------
    collapsed_gene_coverage : Series like
                              Collapsed version
    """
    if '.gz' in gene_coverages:
        gene_coverages_df = pd.read_table(gene_coverages, compression='gzip')
    else:
        gene_coverages_df = pd.read_table(gene_coverages)
    # Not all genes would have the same offset
    # So pad genes with offset < max_offset5p with NAs
    max_offset5p = gene_coverages_df['offset_5p'].max()
    gene_coverages_df[
        'padding_5p_required'] = gene_coverages_df['offset_5p'] - max_offset5p
    zipped_cols = list(
        zip(gene_coverages_df['count'].tolist(), gene_coverages_df[
            'padding_5p_required'].tolist()))
    collapsed = [
        pad_five_prime_or_truncate(
            eval(coverage) / np.nanmean(eval(coverage)), offset_5p,
            target_length) for coverage, offset_5p in zipped_cols
    ]
    metagene = np.nanmean(collapsed, axis=0)
    metagene = pd.Series(
        metagene, index=np.arange(-max_offset5p, len(metagene) - max_offset5p))
    if outfile:
        pickle.dump(metagene, open(outfile, 'wb'))
    else:
        return metagene


def count_reads_per_gene_htseq(bam, bed, prefix=None):
    """Count number of reads following in each region.

    Parameters
    ----------
    bam : str
          Path to bam file (unique mapping only)
    bed : pybedtools.BedTool or str
          Genomic regions to get distance from
    prefix : str
            Prefix to output pickle files

    Returns
    -------
    counts_by_region : Series
                       Series with counts indexed by gene id
    region_lengths : Series
                     Series with gene lengths
    counts_normalized_by_length : Series
                                  Series with normalized counts
    """
    counts_by_region = OrderedDict()
    length_by_region = OrderedDict()
    sorted_bam = HTSeq.BAM_Reader(bam)
    if isinstance(bed, six.string_types):
        bed = pybedtools.BedTool(bed).sort()
    for x, region in enumerate(bed):
        counts = 0
        window = HTSeq.GenomicInterval(
            str(region.chrom).strip(), region.start, region.stop,
            str(region.strand).strip())
        length = region.stop - region.start
        for almnt in sorted_bam[window]:
            counts += 1
        if region.name not in list(counts_by_region.keys()):
            counts_by_region[region.name] = 0
            length_by_region[region.name] = 0
        counts_by_region[region.name] += counts
        length_by_region[region.name] += length
    counts_by_region = pd.Series(counts_by_region)
    length_by_region = pd.Series(length_by_region)
    counts_normalized_by_length = counts_by_region.div(length_by_region)
    if prefix:
        mkdir_p(os.path.dirname(prefix))
        df = pd.concat(
            [counts_by_region, length_by_region, counts_normalized_by_length],
            axis=1)
        df.columns = ['counts', 'length', 'normalized_counts']
        df.to_csv(
            '{}.tsv'.format(prefix), index=True, header=True, sep=str('\t'))
        pickle.dump(counts_by_region,
                    open('{}_counts.pickle'.format(prefix), 'wb'),
                    __PICKLE_PROTOCOL__)
        pickle.dump(length_by_region,
                    open('{}_lengths.pickle'.format(prefix), 'wb'),
                    __PICKLE_PROTOCOL__)
        pickle.dump(counts_normalized_by_length,
                    open('{}_counts_lengths_normalized.pickle'.format(prefix),
                         'wb'), __PICKLE_PROTOCOL__)
    return counts_by_region, length_by_region, counts_normalized_by_length
