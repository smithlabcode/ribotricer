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
from .helpers import list_to_ranges
from .helpers import mkdir_p


def _fix_bed_coltype(bed):
    """Fix bed chrom and name columns to be string

    This is necessary since the chromosome numbers are often interpreted as int
    """
    bed['chrom'] = bed['chrom'].astype(str)
    bed['name'] = bed['name'].astype(str)
    return bed


def collapse_intervals(gene_name, bed, gene_group=None):
    """Extract fasta genewise given coordinates in bed file

    Parameters
    ----------
    gene_name : str
                Gene name
    bed : str
          Path to CDS or 5'UTR or 3'UTR bed
    fasta : str
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

    # Collect all intervals at once
    intervals = zip(gene_group['chrom'], gene_group['start'],
                    gene_group['end'], gene_group['strand'])
    if strand == '+':
        rc = False
    elif strand == '-':
        rc = True

    interval_list = []
    for index, interval in enumerate(intervals):
        strand_interval = interval[3]
        assert strand == strand_interval

        if strand == '+':
            series_range = range(interval[1], interval[2])
        elif strand == '-':
            series_range = range(interval[2], interval[1], -1)

        series = pd.Series(series_range, index=series_range)
        interval_list.append(series)

    interval_combined = interval_list[0]
    for interval_series in interval_list[1:]:
        interval_combined = interval_combined.combine_first(interval_series)
    interval_combined = interval_combined.fillna(0)

    intervals_for_fasta_read = list_to_ranges(interval_combined.index.tolist())
    return intervals_for_fasta_read, strand


def export_gene_fasta(gene_name,
                      bed,
                      chrom_sizes,
                      fasta,
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
    fasta : str
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
    fasta = Fasta(fasta)
    chromsome_lengths = pd.read_table(
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

    # Collect all intervals at once
    intervals = zip(gene_group['chrom'], gene_group['start'],
                    gene_group['end'], gene_group['strand'])
    intervals_list = [list(element) for element in list(intervals)[:]]
    # intervals_list = list(intervals)[:]

    first_interval = intervals_list[0]
    last_interval = intervals_list[-1]
    try:
        chrom_length = chromsome_lengths.loc[str(first_interval[0])]['size']
    except KeyError:
        # for some reason this chromosome is ont part of the bigiwig, so just skip i ([], None)
        raise RuntimeError('chrom length not found')
    # Need to convert to list instead frm tuples
    # TODO fix this?
    # intervals = list(map(list, list(intervals)))
    if strand == '+':
        # For positive strand shift
        # start codon position first_interval[1] by -offset
        if first_interval[1] - offset_5p >= 0:
            first_interval[1] = first_interval[1] - offset_5p
            gene_offset_5p = offset_5p
        else:
            sys.stderr.write('Cannot offset beyond 0 for interval: {}. \
                Set to start of chromsome.\n'.format(first_interval))
            # Reset offset to minimum possible
            gene_offset_5p = first_interval[1]
            first_interval[1] = 0

        if (last_interval[2] + offset_3p <= chrom_length):
            last_interval[2] = last_interval[2] + offset_3p
            gene_offset_3p = offset_3p
        else:
            sys.stderr.write('Cannot offset beyond 0 for interval: {}. \
                             Set to end of chromsome.\n'.format(last_interval))
            gene_offset_3p = chrom_length - last_interval[2]
            # 1-end so chrom_length
            last_interval[2] = chrom_length
    else:
        # Else shift cooridnate of last element in intervals stop by + offset
        if (last_interval[2] + offset_5p <= chrom_length):
            last_interval[2] = last_interval[2] + offset_5p
            gene_offset_5p = offset_5p
        else:
            sys.stderr.write('Cannot offset beyond 0 for interval: {}. \
                             Set to end of chromsome.\n'.format(last_interval))
            gene_offset_5p = chrom_length - last_interval[2]
            # 1-end so chrom_length
            last_interval[2] = chrom_length
        if first_interval[1] - offset_3p >= 0:
            first_interval[1] = first_interval[1] - offset_3p
            gene_offset_3p = offset_3p
        else:
            sys.stderr.write('Cannot offset beyond 0 for interval: {}. \
                Set to start of chromsome.\n'.format(first_interval))
            # Reset offset to minimum possible
            gene_offset_5p = first_interval[1]
            first_interval[1] = 0

    intervals = [tuple(element) for element in intervals_list]
    if strand == '+':
        rc = False
    elif strand == '-':
        rc = True

    interval_list = []
    for index, interval in enumerate(intervals):
        strand_interval = interval[3]
        assert strand == strand_interval

        if strand == '+':
            series_range = range(intervals[index][1], intervals[index][2])
        elif strand == '-':
            series_range = range(intervals[index][2], intervals[index][1], -1)

        series = pd.Series(series_range, index=series_range)
        interval_list.append(series)

    if len(interval_list) == 0:
        # Some genes might not be present in the bigwig at all
        raise RuntimeError('Got empty list! intervals  for chr : {}\n'.format(
            first_interval[0]))

    interval_combined = interval_list[0]
    for interval_series in interval_list[1:]:
        interval_combined = interval_combined.combine_first(interval_series)
    interval_combined = interval_combined.fillna(0)

    intervals_for_fasta_read = list_to_ranges(interval_combined.index.tolist())
    seq = ''
    for interval in intervals_for_fasta_read:
        # Fetching is 1-based for both start and stop and since we are
        # fetching from ranges, there is additional +1
        seq += str(fasta.get_seq(chrom, interval[0] + 1, interval[1] + 1))
    if rc:
        seq = str(Seq(seq, generic_dna).reverse_complement())
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
            gene_offset_5p, gene_offset_3p, seq = export_gene_fasta(
                gene_name, region_bed, chrom_sizes, fasta, gene_group,
                offset_5p, offset_3p)
            with open('{}_{}.fasta'.format(prefix, gene_name),
                      'w') as fh_fasta:
                fh_fasta.write('>{}_5poffset-{}_3poffset-{}\n{}'.format(
                    gene_name, gene_offset_5p, gene_offset_3p, seq))
            outfile.write('{}\t{}\t{}\t{}\n'.format(
                gene_name, int(gene_offset_5p), int(gene_offset_3p), seq))


def get_fasta_sequence(fasta_f, intervals):
    """Extract fasta sequence given a list of intervals.

    Parameters
    ----------
    fasta_f : str
            Path to fasta file

    intervals : list(tuple)
                A list of tuple in the form [(chrom, start, stop, strand)]

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
        # fetching from ranges, there is additional +1
        seq += str(
            fasta.get_seq(chrom, int(interval[1]) + 1, int(interval[2]) + 1))
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
        utr5_coordinates, strand = collapse_intervals(
            gene_name, utr5_bed, utr5_grouped.get_group(gene_name))
        cds_coordinates, strand = collapse_intervals(
            gene_name, cds_bed, cds_grouped.get_group(gene_name))
        utr3_coordinates, strand = collapse_intervals(
            gene_name, utr3_bed, utr3_grouped.get_group(gene_name))

        utr5_intervals = _get_interval(utr5_coordinates, chrom, strand)
        cds_intervals = _get_interval(cds_coordinates, chrom, strand)
        utr3_intervals = _get_interval(utr3_coordinates, chrom, strand)

        utr5seq = get_fasta_sequence(fasta_f, list(utr5_intervals))
        cdsseq = get_fasta_sequence(fasta_f, list(cds_intervals))
        utr3seq = get_fasta_sequence(fasta_f, list(utr3_intervals))
        utr5len = len(utr5seq)
        cdslen = len(cdsseq)
        utr3len = len(utr3seq)

        intervals = list(utr5_intervals) + list(cds_intervals) + list(
            utr3_intervals)
        seq = get_fasta_sequence(fasta_f, intervals)
        with open('{}_{}.fasta'.format(prefix, gene_name), 'w') as fh_fasta:
            fh_fasta.write('>{}_utr5len={};cdslen={};utr3len={}\n{}'.format(
                gene_name, utr5len, cdslen, utr3len, seq))
