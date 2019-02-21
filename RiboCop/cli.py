"""Command line interface for RiboCop
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import click
import sys
import os
from click_help_colors import HelpColorsGroup
from . import __version__
from .detect_orfs import detect_orfs
from .prepare_orfs import prepare_orfs
from .count_orfs import count_orfs

click.disable_unicode_literals_warning = True
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(
    cls=HelpColorsGroup,
    help_headers_color='yellow',
    help_options_color='green')
@click.version_option(version=__version__)
def cli():
    """RiboCop: Tool for detecting translating ORF from Ribo-seq data"""
    pass


###################### prepare-orfs function #########################################
@cli.command(
    'prepare-orfs',
    context_settings=CONTEXT_SETTINGS,
    help='Extract candidate ORFS based on GTF and FASTA files')
@click.option('--gtf', help='Path to GTF file', required=True)
@click.option('--fasta', help='Path to FASTA file', required=True)
@click.option('--prefix', help='Prefix to output file', required=True)
@click.option(
    '--min_orf_length',
    type=int,
    default=60,
    show_default=True,
    help='The minimum length (nts) of ORF to include')
@click.option(
    '--start_codons',
    default='ATG',
    show_default=True,
    help='Comma separated list of start codons')
@click.option(
    '--stop_codons',
    default='TAG,TAA,TGA',
    show_default=True,
    help='Comma separated list of stop codons')
@click.option(
    '--longest',
    help='Choose the most upstream start codon if multiple in frame ones exist',
    is_flag=True)
def prepare_orfs_cmd(gtf, fasta, prefix, min_orf_length, start_codons,
                     stop_codons, longest):
    if not os.path.isfile(gtf):
        sys.exit('Error: GTF file is not found')

    if not os.path.isfile(fasta):
        sys.exit('Error: FASTA file is not found')

    if min_orf_length <= 0:
        sys.exit('Error: min ORF length at least to be 1')

    start_codons = set(
        [x.strip().upper() for x in start_codons.strip().split(',')])
    if not start_codons:
        sys.exit('Error: start codons cannot be empty')
    if not all(
        [len(x) == 3 and set(x) <= {'A', 'C', 'G', 'T'}
         for x in start_codons]):
        sys.exit('Error: invalid codon, only A, C, G, T allowed')

    stop_codons = set(
        [x.strip().upper() for x in stop_codons.strip().split(',')])
    if not stop_codons:
        sys.exit('Error: stop codons cannot be empty')
    if not all(
        [len(x) == 3 and set(x) <= {'A', 'C', 'G', 'T'} for x in stop_codons]):
        sys.exit('Error: invalid codon, only A, C, G, T allowed')

    prepare_orfs(gtf, fasta, prefix, min_orf_length, start_codons, stop_codons,
                 longest)


###################### detect-orfs function #########################################
@cli.command(
    'detect-orfs',
    context_settings=CONTEXT_SETTINGS,
    help='Detect translating ORFs from BAM file')
@click.option('--bam', help='Path to BAM file', required=True)
@click.option(
    '--ribocop_index',
    help=('Path to the index file of RiboCop\n'
          'This file should be generated using RiboCop prepare-orfs'),
    required=True)
@click.option('--prefix', help='Prefix to output file', required=True)
@click.option(
    '--stranded',
    type=click.Choice(['yes', 'no', 'reverse']),
    default=None,
    show_default=True,
    help=
    ('whether the data is from a strand-specific assay'
     ' If not provided, the experimental protocol will be automatically inferred'
     ))
@click.option(
    '--read_lengths',
    default=None,
    show_default=True,
    help=('Comma separated read lengths to be used, such as 28,29,30\n'
          'If not provided, it will be automatically determined by assessing'
          ' the metagene periodicity'))
@click.option(
    '--psite_offsets',
    default=None,
    show_default=True,
    help=('Comma separated P-site offsets for each read length '
          'matching the read lengths provided.\n'
          'If not provided, reads from different read lengths will be '
          'automatically aligned using cross-correlation'))
@click.option(
    '--report_all',
    help=('Whether output all ORFs including those '
          'non-translating ones'),
    is_flag=True)
def detect_orfs_cmd(bam, ribocop_index, prefix, stranded, read_lengths,
                    psite_offsets, report_all):
    if not os.path.isfile(bam):
        sys.exit('Error: BAM file is not found')

    if not os.path.isfile(ribocop_index):
        sys.exit('Error: RiboCop index file is not found')

    if read_lengths is not None:
        try:
            read_lengths = [
                int(x.strip()) for x in read_lengths.strip().split(',')
            ]
        except:
            sys.exit('Error: cannot convert read_lengths into integers')
        if not all([x > 0 for x in read_lengths]):
            sys.exit('Error: read length must be positive')

    if read_lengths is None and psite_offsets is not None:
        sys.exit(
            'Error: psite_offsets only allowed when read_lengths is provided')
    if read_lengths is not None and psite_offsets is not None:
        try:
            psite_offsets = [
                int(x.strip()) for x in psite_offsets.strip().split(',')
            ]
        except:
            sys.exit('Error: cannot convert psite_offsets into integers')
        if len(read_lengths) != len(psite_offsets):
            sys.exit('Error: psite_offsets must match read_lengths')
        if not all(x >= 0 for x in psite_offsets):
            sys.exit('Error: P-site offset must be >= 0')
        if not all(x > y for (x, y) in zip(read_lengths, psite_offsets)):
            sys.exit('Error: P-site offset must be smaller than read length')
        psite_offsets = dict(zip(read_lengths, psite_offsets))
    if stranded == 'yes':
        stranded = 'forward'
    detect_orfs(bam, ribocop_index, prefix, stranded, read_lengths,
                psite_offsets, report_all)


###################### count-orfs function #########################################
@cli.command(
    'count-orfs',
    context_settings=CONTEXT_SETTINGS,
    help='Count reads for detected ORFs')
@click.option(
    '--ribocop_index',
    help=('Path to the index file of RiboCop\n'
          'This file should be generated using RiboCop prepare-orfs'),
    required=True)
@click.option(
    '--detected_orfs',
    help=('Path to the detected orfs file\n'
          'This file should be generated using RiboCop detect-orfs'),
    required=True)
@click.option('--features', help='set of ORF types')
@click.option('--prefix', help='Prefix to output file', required=True)
@click.option(
    '--report_all',
    help=('Whether output all ORFs including those '
          'non-translating ones'),
    is_flag=True)
def detect_orfs_cmd(bam, ribocop_index, prefix, stranded, read_lengths,
                    psite_offsets, report_all):
    if not os.path.isfile(bam):
        sys.exit('Error: BAM file is not found')

    if not os.path.isfile(ribocop_index):
        sys.exit('Error: RiboCop index file is not found')

    if read_lengths is not None:
        try:
            read_lengths = [
                int(x.strip()) for x in read_lengths.strip().split(',')
            ]
        except:
            sys.exit('Error: cannot convert read_lengths into integers')
        if not all([x > 0 for x in read_lengths]):
            sys.exit('Error: read length must be positive')

    if read_lengths is None and psite_offsets is not None:
        sys.exit(
            'Error: psite_offsets only allowed when read_lengths is provided')
    if read_lengths is not None and psite_offsets is not None:
        try:
            psite_offsets = [
                int(x.strip()) for x in psite_offsets.strip().split(',')
            ]
        except:
            sys.exit('Error: cannot convert psite_offsets into integers')
        if len(read_lengths) != len(psite_offsets):
            sys.exit('Error: psite_offsets must match read_lengths')
        if not all(x >= 0 for x in psite_offsets):
            sys.exit('Error: P-site offset must be >= 0')
        if not all(x > y for (x, y) in zip(read_lengths, psite_offsets)):
            sys.exit('Error: P-site offset must be smaller than read length')
        psite_offsets = dict(zip(read_lengths, psite_offsets))
    if stranded == 'yes':
        stranded = 'forward'
    detect_orfs(bam, ribocop_index, prefix, stranded, read_lengths,
                psite_offsets, report_all)
