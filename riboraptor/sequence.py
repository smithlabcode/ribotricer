"""Utilities for extracting sequence from fasta.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import warnings
import os
import re
import sys
import pybedtools
import pandas as pd
import numpy as np

from .genome import _get_sizes
from .genome import _get_bed
from .genome import __GENOMES_DB__

from .helpers import merge_intervals

from .interval import Interval
from .fasta import FastaReader


def gene_sequence(gene_group, fasta, offset_5p=0, offset_3p=0):
    """Extract seq genewise given coordinates in bed file

    Parameters
    ----------
    gene_group: DataFrame
                gene group from bed file
    fasta  str
           Path to fasta file
    offset_5p: int (positive)
               Number of bases to count upstream (5')
    offset_3p: int (positive)
               Number of bases to count downstream (3')

    Returns
    -------
    sequence: str
              sequence of the gene
    gene_offset_5p: Gene wise 5 prime offset
                    This might be different from `offset_5p` in cases where
                    `offset_5p` leads to a negative coordinate
    gene_offset_3p: Gene wise 3 prime offset
                    This might be different from `offset_3p` in cases where
                    `offset_3p` leads to position beyond chromsome length
    """
    if offset_5p < 0 or offset_3p < 0:
        raise RuntimeError('Offsets must be non-negative')
        sys.exit(1)
    if not isinstance(fasta, FastaReader):
        fasta = FastaReader(fasta)
    chromosome_lengths = fasta.chromosomes

    if len(gene_group['strand'].unique()) != 1:
        raise RuntimeError('Multiple strands?: {}'.format(gene_group))
    if len(gene_group['chrom'].unique()) != 1:
        raise RuntimeError('Chromosome not unique for: {}'.format(gene_group))
    chrom = gene_group['chrom'].unique()[0]
    strand = gene_group['strand'].unique()[0]

    # convert zero-based half-open interval to one-based full-closed
    intervals = list(
        zip(gene_group['chrom'], gene_group['start'] + 1, gene_group['end'],
            gene_group['strand']))

    intervals = [Interval(i[0], i[1], i[2], i[3]) for i in intervals]
    intervals_combined, gene_offset_5p, gene_offset_3p = merge_intervals(
        intervals, chromosome_lengths, offset_5p, offset_3p, False)
    sequences = fasta.query(intervals_combined)
    sequence = ''.join(sequences)
    if strand == '-':
        sequence = fasta.reverse_complement(sequence)
    return (sequence, gene_offset_5p, gene_offset_3p)


def export_gene_sequences(bed, fasta, saveto=None, offset_5p=0, offset_3p=0):
    """Export all gene sequences.

    Parameters
    ----------
    bed: str
         Path to CDS or 5'UTR or 3'UTR bed
    fasta: str
           Path to fasta to fetch the sequences from
    saveto: str
            Path to write output tsv file
    offset_5p: int (positive)
               Number of bases to count upstream (5')
    offset_3p: int (positive)
               Number of bases to count downstream (3')

    Returns
    -------
    gene_profiles: file
                   with the following format:
                   gene1\t5poffset1\t3poffset1\tsequence\n
                   gene2\t5poffset2\t3poffset2\tsequence\n
    """
    if bed.lower().split('_')[0] in __GENOMES_DB__:
        genome, region_type = bed.lower().split('_')
        bed = _get_bed(region_type, genome)
    bed_df = pybedtools.BedTool(bed).sort().to_dataframe()
    bed_df['chrom'] = bed_df['chrom'].astype(str)
    bed_df['name'] = bed_df['name'].astype(str)
    bed_grouped = bed_df.groupby('name')
    to_write = 'gene_name\toffset_5p\toffset_3p\tsequence\n'
    for gene_name, gene_group in bed_grouped:
        sequence, gene_offset_5p, gene_offset_3p = gene_sequence(
            gene_group, fasta, offset_5p, offset_3p)
        to_write += '{}\t{}\t{}\t{}\n'.format(gene_name,
                                              int(gene_offset_5p),
                                              int(gene_offset_3p),
                                              sequence)
    if not saveto:
        return to_write
    else:
        with open(saveto, 'w') as outfile:
            outfile.write(to_write)
