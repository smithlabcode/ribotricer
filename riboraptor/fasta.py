from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os
import re
import sys
import pybedtools
import pandas as pd
import numpy as np

from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from .genome import _get_sizes
from .genome import _get_bed
from .genome import __GENOMES_DB__

from .helpers import collapse_bed_intervals
from .helpers import list_to_ranges
from .helpers import mkdir_p
from .helpers import _fix_bed_coltype


def get_fasta_sequence(fasta_f, intervals):
    """Extract fasta sequence given a list of intervals.

    Parameters
    ----------
    fasta_f : str
            Path to fasta file

    intervals : list(tuple)
                A list of tuple in the form [(chrom, start, stop, strand)]
                NOTE: 1-based start and stop only!

    Returns
    -------
    seq : list
          List of sequences at intervals
    """
    fasta = Fasta(fasta_f)
    chrom = intervals[0][0]
    strand = intervals[0][-1]
    seq = ''
    for interval in intervals:
        # Fetching is 1-based for both start and stop and since we are
        # fetching from ranges,
        # This is already accounted for internally in the pipeline
        # when using the method ""there is additional +1
        seq += str(fasta.get_seq(chrom, int(interval[1]), int(interval[2])))
    if strand == '-':
        seq = str(Seq(seq, generic_dna).reverse_complement())
    return seq


def _get_interval(coordinates, chrom, strand):
    intervals = []
    for coord in coordinates:
        start = coord[0]
        end = coord[1]
        intervals.append((chrom, start, end, strand))
    return intervals


def export_fasta_from_bed(gene_name,
                          bed,
                          chrom_sizes,
                          fasta_f,
                          gene_group=None,
                          offset_5p=0,
                          offset_3p=0):
    """Extract fasta genewise given coordinates in bed file

    Parameters
    ----------
    gene_name : str
                Gene name
    bed : str
          Path to CDS or 5'UTR or 3'UTR bed
    fasta_f : str
            Path to fasta file
    chrom_sizes : str
                  Path to chrom.sizes file
    offset_5p : int (positive)
                Number of bases to count upstream (5')
    offset_3p : int (positive)
                Number of bases to count downstream (3')

    Returns
    -------
    gene_offset : int
                Gene wise offsets
    """
    assert (offset_5p >= 0)
    assert (offset_3p >= 0)
    chromosome_lengths = pd.read_table(
        chrom_sizes, names=['chrom', 'size']).set_index('chrom')
    if not isinstance(bed, pd.DataFrame):
        bed = pybedtools.BedTool(bed).to_dataframe()
        bed = _fix_bed_coltype(bed)
    assert gene_name in bed['name'].tolist()
    if gene_group is None:
        gene_group = bed[bed['name'] == gene_name]

    assert len(gene_group['strand'].unique()) == 1
    assert len(gene_group['chrom'].unique()) == 1
    chrom = gene_group['chrom'].unique()[0]
    strand = gene_group['strand'].unique()[0]
    if strand == '+':
        rc = False
    elif strand == '-':
        rc = True

    # Collect all intervals at once
    intervals = list(
        zip(gene_group['chrom'], gene_group['start'], gene_group['end'],
            gene_group['strand']))
    query_intervals, fasta_onebased_intervals, gene_offset_5p, gene_offset_3p = collapse_bed_intervals(
        intervals, chromosome_lengths, offset_5p, offset_3p)

    seq = get_fasta_sequence(fasta_f, list(fasta_onebased_intervals))
    return gene_offset_5p, gene_offset_3p, str(seq)


def export_all_fasta(region_bed_f,
                     chrom_sizes,
                     fasta,
                     prefix,
                     offset_5p=60,
                     offset_3p=0,
                     ignore_tx_version=True):
    """Export all gene coverages.

    Parameters
    ----------
    region_bed_f : str
                   Path to region bed file (CDS/3'UTR/5'UTR)
                   with bed name column as gene
    chrom_sizes : str
                  Path to chrom.sizes file
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
    """
    if region_bed_f.lower().split('_')[0] in __GENOMES_DB__:

        genome, region_type = region_bed_f.lower().split('_')
        region_bed_f = _get_bed(region_type, genome)

    region_bed = pybedtools.BedTool(region_bed_f).sort().to_dataframe()
    region_bed = _fix_bed_coltype(region_bed)
    # Group intervals by gene name
    cds_grouped = region_bed.groupby('name')

    mkdir_p(os.path.dirname(prefix))
    with open('{}_all_genes.tsv'.format(prefix), 'w') as outfile:
        outfile.write('gene_name\toffset_5p\toffset_3p\tseq\n')
        for gene_name, gene_group in cds_grouped:
            if ignore_tx_version:
                gene_name = re.sub(r'\.[0-9]+', '', gene_name)
            gene_offset_5p, gene_offset_3p, seq = export_fasta_from_bed(
                gene_name, region_bed, chrom_sizes, fasta, gene_group,
                offset_5p, offset_3p)
            with open('{}_{}.fasta'.format(prefix, gene_name),
                      'w') as fh_fasta:
                fh_fasta.write('>{}_5poffset-{}_3poffset-{}\n{}'.format(
                    gene_name, gene_offset_5p, gene_offset_3p, seq))
            outfile.write('{}\t{}\t{}\t{}\n'.format(
                gene_name, int(gene_offset_5p), int(gene_offset_3p), seq))


def complete_gene_fasta(utr5_bed_f, cds_bed_f, utr3_bed_f, fasta_f, prefix):
    """Merge Utr5, CDS, UTR3 coordinates to get one fasta.

    Parameters
    ----------
    utr5_bed : str
               Path to 5'UTR bed
    cds_bed : str
              Path to CDS bed
    utr3_bed : str
               Path to 3'UTR bed
    """
    utr5_bed = _fix_bed_coltype(
        pybedtools.BedTool(utr5_bed_f).sort().to_dataframe())
    cds_bed = _fix_bed_coltype(
        pybedtools.BedTool(cds_bed_f).sort().to_dataframe())
    utr3_bed = _fix_bed_coltype(
        pybedtools.BedTool(utr3_bed_f).sort().to_dataframe())
    # Group intervals by gene namei

    utr5_grouped = utr5_bed.groupby('name')
    cds_grouped = cds_bed.groupby('name')
    utr3_grouped = utr3_bed.groupby('name')

    mkdir_p(os.path.dirname(prefix))

    # For now just write records where wll three annotations are complete
    utr5_genes = set(list(utr5_grouped.groups))
    cds_genes = set(list(cds_grouped.groups))
    utr3_genes = set(list(utr3_grouped.groups))

    common_genes = utr5_genes.intersection(cds_genes).intersection(utr3_genes)

    for gene_name in common_genes:
        chrom = cds_grouped.get_group(gene_name)['chrom'].unique()[0]
        strand = cds_grouped.get_group(gene_name)['strand'].unique()[0]

        utr5_group = utr5_grouped.get_group(gene_name)
        cds_group = cds_grouped.get_group(gene_name)
        utr3_group = utr3_grouped.get_group(gene_name)

        utr5_intervals = list(
            zip(utr5_group['chrom'], utr5_group['start'], utr5_group['end'],
                utr5_group['strand']))
        cds_intervals = list(
            zip(cds_group['chrom'], cds_group['start'], cds_group['end'],
                cds_group['strand']))
        utr3_intervals = list(
            zip(utr3_group['chrom'], utr3_group['start'], utr3_group['end'],
                utr3_group['strand']))

        _, utr5_onebased_intervals, _, _ = collapse_bed_intervals(
            utr5_intervals, chromosome_lengths, offset_5p, offset_3p)
        _, cds_onebased_intervals, _, _ = collapse_bed_intervals(
            cds_intervals, chromosome_lengths, offset_5p, offset_3p)
        _, utr3_onebased_intervals, _, _ = collapse_bed_intervals(
            utr3_intervals, chromosome_lengths, offset_5p, offset_3p)

        utr5seq = get_fasta_sequence(fasta_f, list(utr5_onebased_intervals))
        cdsseq = get_fasta_sequence(fasta_f, list(cds_onebased_intervals))
        utr3seq = get_fasta_sequence(fasta_f, list(utr3_onebased_intervals))
        utr5len = len(utr5seq)
        cdslen = len(cdsseq)
        utr3len = len(utr3seq)

        intervals = list(utr5_onebased_intervals) + list(
            cds_onebased_intervals) + list(utr3_onebased_intervals)
        seq = get_fasta_sequence(fasta_f, intervals)
        with open('{}_{}.fasta'.format(prefix, gene_name), 'w') as fh_fasta:
            fh_fasta.write('>{}_utr5len={};cdslen={};utr3len={}\n{}'.format(
                gene_name, utr5len, cdslen, utr3len, seq))
