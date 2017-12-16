"""Utilities for read counting operations."""
from __future__ import division
from collections import Counter
from collections import defaultdict

import HTSeq
import pandas as pd
import pysam
import pybedtools
from pyfaidx import Fasta

import numpy as np
import re
import sys

from.wig import WigReader


def mapped_reads_count(bam):
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
    if isinstance(bam, str):
        bam = pysam.AlignmentFile(bam, 'rb')
    n_mapped = len([query for query in bam.fetch()
                    if query.mapping_quality == 255
                    and not query.is_secondary])
    return n_mapped

def fragments_length_distribution(bam):
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
    if isinstance(bam, str):
        bam = pysam.AlignmentFile(bam, 'rb')
    return Counter([r.query_length for r in bam.fetch()])

def get_coverage(bam):
    if isinstance(bam, str):
        bam = pysam.AlignmentFile(bam, 'rb')

    coverage = defaultdict(Counter)
    for read in bam.fetch():
        if read.is_secondary:
            continue
        if read.is_reverse:
            strand = '-'
        else:
            strand = '+'
        reference_pos = read.get_reference_positions()
        if strand == '+':
            position = reference_pos[0]  #[offsets[read_length]-1]
        else:
            ## Negative strand so no need to adjust for 0-based indexing
            position = reference_pos[-1]  #[-offsets[read_length]]
        query_length = read.query_length
        coverage['{}:{}'.format(read.reference_name, position)][query_length] += 1
    return coverage


def get_gene_coverage(gene_name, bed, bw, master_offset=0):
    """Get gene coverage

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
    #if isinstance(bw, str):
    bw = WigReader(bw)

    chromsome_lengths = bw.get_chromosomes

    bed = pybedtools.BedTool(bed).to_dataframe()
    assert gene_name in bed['name'].tolist()
    gene_group = bed[bed['name']==gene_name]

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
        sys.stderr.write('Got empty list! intervals  for chr : {}\n'.format(intervals[0][0]))
        return ([], None)

    coverage_combined = interval_coverage_list[0]
    for interval_coverage in interval_coverage_list[1:]:
        coverage_combined = coverage_combined.combine_first(interval_coverage)
    coverage_combined = coverage_combined.fillna(0)
    index_to_genomic_pos_map = pd.Series(coverage_combined.index.tolist(), index=np.arange(len(coverage_combined))-gene_offset)
    intervals_for_fasta_query = []
    for pos in index_to_genomic_pos_map.values:
        intervals_for_fasta_query.append((chrom, pos, pos+1, strand))
    coverage_combined = coverage_combined.reset_index(drop=True)
    coverage_combined = coverage_combined.rename(lambda x: x - gene_offset)

    return (coverage_combined, intervals_for_fasta_query,
            index_to_genomic_pos_map, gene_offset)

