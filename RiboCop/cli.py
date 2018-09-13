"""Command line interface for riboraptor
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import sys
from textwrap import dedent

import click
from click_help_colors import HelpColorsGroup
from . import __version__
from .infer_protocol import infer_protocol
from .orf import prepare_orfs

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
    help='extract putative orfs based on GTF and FASTA files')
@click.option('--gtf', help='Path to GTF file')
@click.option('--fasta', help='Path to FASTA file')
@click.option('--prefix', help='Prefix to output file')
def prepare_orfs_cmd(gtf, fasta, prefix):
    prepare_orfs(gtf, fasta, prefix)


###################### infer-protocol function #########################################
@cli.command(
    'infer-protocol',
    context_settings=CONTEXT_SETTINGS,
    help='Infer protocol from BAM')
@click.option('--bam', help='Path to bam file')
@click.option(
    '--gtf', help='Path to gtf file to be used for defining the gene strands')
@click.option('--prefix', help='Prefix to output protocol file')
@click.option(
    '--n_reads',
    type=int,
    default=20000,
    help='Number of mapped reads to use for estimation')
def infer_protocol_cmd(bam, gtf, prefix, n_reads):
    protocol, forward_mapped, reverse_mapped = infer_protocol(
        bam, gtf, prefix, n_reads)
    total = forward_mapped + reverse_mapped
    print(
        dedent('''\
                 Forward mapped proportion: {} ({:.4f})
                 Reverse mapped proportion: {} ({:.4f})
                 Likely protocol: {}'''
               .format(forward_mapped, forward_mapped / total, reverse_mapped,
                       reverse_mapped / total, protocol)))
