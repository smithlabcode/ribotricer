"""Utilities for read counting operations."""
from __future__ import division
from collections import Counter
from collections import defaultdict

import HTSeq
import pandas as pd
import pysam
import pybedtools
import pyBigWig
from pyfaidx import Fasta

import numpy as np
import warnings
import re
import sys

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

"""
def collapse_region(records):
    record = records[0]
    chrom = record['chrom']
    strand = record['strand']
    region = pd.Series(range(record['start'], record['end']))
    for record in records:
        region = region.combine_first(pd.Series(range(record['start'], record['end'])))
    return (chrom, region[0], region[-1], strand)

def get_coordinates(bed, gene):
    '''Get collapsed coordiantes of a gene from the given bed file'''
    bed = pybedtools.BedTool(bed).to_dataframe()
    assert gene in bed['name']
    region = collapse_region(bed[bed['name']==gene].T.to_dict().values())
    return region
"""

class WigReader(object):
    """Class for reading and querying wigfiles"""

    def __init__(self, wig_location):
        """
        Arguments
        ---------
        wig_location: Path to bigwig
        """
        self.wig_location = wig_location
        try:
            self.wig = pyBigWig.open(self.wig_location)
        except Exception as e:
            raise Exception('Error reading wig file: {}'.format(e))

    def query(self, intervals):
        """ Query regions for scores
        Arguments
        ---------
        intervals: list of tuples
            A list of tuples with the following format: (chr, chrStart, chrEnd, strand)
        Returns
        -------
        scores: np.array
            A numpy array containing scores for each tuple
        """
        scores = []
        chrom_lengths = self.get_chromosomes
        for chrom, chromStart, chromEnd, strand in intervals:
            if chrom not in list(chrom_lengths.keys()):
                warnings.warn(
                    'Chromosome {} does not appear in the bigwig'.format(chrom), UserWarning)
                continue

            chrom_length = chrom_lengths[chrom]
            if int(chromStart) > chrom_length:
                raise Exception('Chromsome start point exceeds chromosome length: {}>{}'.format(
                    chromStart, chrom_length))
            elif int(chromEnd) > chrom_length:
                raise Exception('Chromsome end point exceeds chromosome length: {}>{}'.format(
                    chromEnd, chrom_length))
            score = self.wig.values(chrom, int(chromStart), int(chromEnd))
            if strand == '-':
                score.reverse()
            scores.append(score)
        return np.array(scores)

    @property
    def get_chromosomes(self):
        """Return list of chromsome and their sizes
        as in the wig file
        Returns
        -------
        chroms: dict
            Dictionary with {"chr": "Length"} format
        """
        return self.wig.chroms()

def get_gene_coverage(gene_name, bed, bw, master_offset=0):
    """Get gene coverage

    Parameters
    ----------

    gene_name: str
        Gene name
    bed: str
        Path to CDS or 5'UTR or 3'UTR bed
    bw: str
        Path to bigwig to fetch the scores from
    master_offset: int
        How much bases upstream?

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

def get_fasta_sequence(fasta, intervals):
    if isinstance(fasta, str):
        fasta = Fasta(fasta)
    sequence = []
    for interval in intervals:
        chrom, start, stop, strand = interval
        if strand == '+':
            seq = fasta[chrom][int(start):int(stop)].seq
        elif strand == '-':
            seq = fasta[chrom][int(start):int(stop)].reverse.seq
        sequence.append(seq)
    return ('').join(sequence)

def get_region_sizes(region_bed_grouped):
    region_sizes = {}
    for gene_name, gene_group in region_bed_grouped:
        ## Get rid of trailing dots
        gene_name = re.sub(r'\.[0-9]+', '', gene_name)
        # Collect all intervals at once
        intervals = zip(gene_group['chrom'], gene_group['start'],
                        gene_group['end'], gene_group['strand'])
        for interval in intervals:
            if gene_name not in region_sizes:
                # End is always 1-based so does not require +1
                region_sizes[gene_name] = interval[2] - interval[1]
            else:
                region_sizes[gene_name] += interval[2] - interval[1]
    return region_sizes

def sort_genes_lengthwise(bed):
    bed = pybedtools.BedTool(bed).to_dataframe()
    bed_grouped = bed.groupby('name')
    sizes = get_region_sizes(bed_grouped)
    sizes_df = pd.DataFrame(sizes.items(), columns=['name', 'length'])
    
    sizes_df.sort_values(by='length', ascending=False, inplace=True)
    return [tuple(x) for x in sizes_df[['name', 'length']].values]
