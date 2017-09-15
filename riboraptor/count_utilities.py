"""Utilities for read counting operations."""
from __future__ import division
from collections import Counter
from collections import defaultdict

import HTSeq
import pandas as pd
import pysam
import pybedtools

def bed_to_genomic_interval(bed):
    '''
        Converts bed file to genomic interval (htseq format) file
    '''
    for interval in bed:
        yield HTSeq.GenomicPosition(
            str(interval.chrom), interval.start, str(interval.strand)
        )


def count_mapped(bam):
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

def fragments_lengths(bam):
    """Count fragment lengths.

    Parameters
    ----------

    bam : str
        Path

    Returns
    -------
    lengths : array_like
    """
    if isinstance(bam, str):
        bam = pysam.AlignmentFile(bam, 'rb')
    return [r.query_length for r in bam.fetch()]

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

def get_closest(bam, regions, half_window_width=100):
    '''
    Given a bam file and a list of regions returns a dataframe with the distance of each read from the closest region
    bam -- bam file
    regions -- bed file, genomic regions to get distance from
    half_window_width -- int, distance around region to record
    '''
    profile = {}
    sorted_bam = HTSeq.BAM_Reader(bam.fn)
    for x, tss in enumerate(bed_to_genomic_interval(regions)):
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


def start_stop_codon_reads(bam, start_codon_bed, stop_codon_bed):
    """Return reads around start and stop codons """
    start_closest = get_closest(pybedtools.BedTool(bam), start_codon_bed)
    stop_closest = get_closest(pybedtools.BedTool(bam), stop_codon_bed)
    return start_closest, stop_closest

def get_coverage_near_start_codons(start_codon_bed,
                                   coverage_dict,
                                   n_nucleotides=600):

    bed_df = start_codon_bed.to_dataframe()
    pos_names = bed_df.chrom.str.cat(bed_df.start.astype(str), sep=':').tolist()

    coverage = defaultdict(Counter)

    for start_codon in pos_names:
        chrom, start_position = start_codon.split(':')
        start_position = int(start_position)
        cur_pointer = -n_nucleotides
        while cur_pointer < n_nucleotides:
            curr_pos = start_position + cur_pointer
            counts_counter = coverage_dict['{}:{}'.format(chrom, curr_pos)]
            coverage[cur_pointer] += counts_counter
            cur_pointer += 1
    return pd.Series(coverage)


def summarize_counts(coverage, length_range):
    """Summarize counts over a length range"""
    coverage_series = coverage.copy()
    for index, values in coverage_series.iteritems():
        coverage_series[index] = pd.Series(values)[length_range].sum()
    return coverage_series
