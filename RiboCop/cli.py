"""Command line interface for riboraptor
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import sys
from textwrap import dedent
import numpy as np
import pandas as pd

import click
from click_help_colors import HelpColorsGroup
from . import __version__
from .infer_protocol import infer_protocol
from .orf import *
from .utils import parse_ccds

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
    help='Extract putative orfs based on GTF and FASTA files')
@click.option('--gtf', help='Path to GTF file')
@click.option('--fasta', help='Path to FASTA file')
@click.option('--prefix', help='Prefix to output file')
def prepare_orfs_cmd(gtf, fasta, prefix):
    prepare_orfs(gtf, fasta, prefix)


###################### detect-orfs function #########################################
@cli.command(
    'detect-orfs',
    context_settings=CONTEXT_SETTINGS,
    help='Detect translating ORFs from BAM file')
@click.option('--bam', help='Path to BAM file')
@click.option('--prefix', help='Prefix to output file')
@click.option('--gtf', help='Path to GTF file')
@click.option('--annotation', help='Path to annotation file')
def detect_orfs_cmd(bam, prefix, gtf, annotation):
    detect_orfs(bam, prefix, gtf=gtf, annotation=annotation)


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
    infer_protocol(bam, gtf, prefix, n_reads)


###################### parse-ccds function #########################################
@cli.command(
    'parse-ccds',
    context_settings=CONTEXT_SETTINGS,
    help='Generate periodicity for ccds')
@click.option('--annotation', help='Path to annotation file')
@click.option('--orfs', help='Path to ORF file')
@click.option('--saveto', help='Path of output file')
def parse_ccds_cmd(annotation, orfs, saveto):
    parse_ccds(annotation, orfs, saveto)

###################### parse-annotation function #########################################
@cli.command(
    'parse-annotation',
    context_settings=CONTEXT_SETTINGS,
    help='Parse annotation file to extract putative ORFs')
@click.option('--annotation', help='Path to annotation file')
def parse_annotation_cmd(annotation):
    parse_annotation(annotation)


###################### plot-read-lengths function #########################################
@cli.command(
    'plot-metagene',
    context_settings=CONTEXT_SETTINGS,
    help='Plot read length distribution')
@click.option('--prefix', help='Prefix to output file')
def plot_metagene_cmd(prefix):
    metagenes = {}
    s = pd.Series(np.random.randint(30, size=500), index=np.arange(-50, 450))
    read_lengths = {28: 100023434, 27: 2690000, 30: 1234908, 31: 790000, 32: 8000, 26: 123490}
    for length in read_lengths:
        metagenes[length] = s
    plot_metagene(metagenes, read_lengths, prefix, 100)


