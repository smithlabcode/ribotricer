"""Command line interface for RiboCop
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import pandas as pd

import click
from click_help_colors import HelpColorsGroup
from . import __version__
from .detect_orfs import detect_orfs
from .infer_protocol import infer_protocol
from .prepare_orfs import prepare_orfs
from .utils import parse_ccds
from .utils import test_periodicity
from .utils import benchmark
from .utils import theta_dist
from .utils import theta_rna

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
    help='Extract candidate orfs based on GTF and FASTA files')
@click.option('--gtf', help='Path to GTF file', required=True)
@click.option('--fasta', help='Path to FASTA file', required=True)
@click.option('--prefix', help='Prefix to output file', required=True)
@click.option(
    '--region',
    help='Regions to extract, comma separated str (e.g.  "cds,utr")')
def prepare_orfs_cmd(gtf, fasta, prefix, region):
    region = [x.strip() for x in region.strip().split(',')]
    prepare_orfs(gtf, fasta, prefix, region)


###################### detect-orfs function #########################################
@cli.command(
    'detect-orfs',
    context_settings=CONTEXT_SETTINGS,
    help='Detect translating ORFs from BAM file')
@click.option('--bam', help='Path to BAM file', required=True)
@click.option('--prefix', help='Prefix to output file', required=True)
@click.option('--gtf', help='Path to GTF file', required=True)
@click.option('--annotation', help='Path to annotation file', required=True)
def detect_orfs_cmd(bam, prefix, gtf, annotation):
    detect_orfs(bam, prefix, gtf=gtf, annotation=annotation)


###################### export-rna function #########################################
@cli.command(
    'export-rna',
    context_settings=CONTEXT_SETTINGS,
    help='Detect translating ORFs from RNA BAM file')
@click.option('--bam', help='Path to BAM file')
@click.option('--prefix', help='Prefix to output file')
@click.option('--gtf', help='Path to GTF file')
@click.option('--annotation', help='Path to annotation file')
def detect_orfs_cmd(bam, prefix, gtf, annotation):
    detect_orfs(bam, prefix, gtf=gtf, annotation=annotation, testRNA=True)


###################### infer-protocol function #########################################
@cli.command(
    'infer-protocol',
    context_settings=CONTEXT_SETTINGS,
    help='Infer protocol from BAM')
@click.option('--bam', help='Path to bam file', required=True)
@click.option(
    '--gtf', help='Path to gtf file to be used for defining the gene strands', required=True)
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


###################### test-periodicity function #########################################
@cli.command(
    'test-periodicity',
    context_settings=CONTEXT_SETTINGS,
    help='Test different method for periodicity score')
@click.option('--orf', help='Path to ORF file')
@click.option('--prefix', help='Prefix of output')
@click.option('--method', help='Method for periodicity score')
def test_periodicity_cmd(orf, prefix, method):
    test_periodicity(orf, prefix, method)


###################### benchmark function #########################################
@cli.command(
    'benchmark',
    context_settings=CONTEXT_SETTINGS,
    help='Test different method for periodicity score')
@click.option('--rna', help='Path to rna ORF file', required=True)
@click.option('--ribo', help='Path to ribo ORF file', required=True)
@click.option('--prefix', help='Prefix to output file', required=True)
@click.option(
    '--cutoff', type=int, default=5, help='Cutoff of number of reads')
def benchmark_cmd(rna, ribo, prefix, cutoff):
    benchmark(rna, ribo, prefix, cutoff)


###################### theta distribution function #########################################
@cli.command(
    'theta-dist',
    context_settings=CONTEXT_SETTINGS,
    help='Test distribution of angles')
@click.option('--rna', help='Path to rna ORF file')
@click.option('--ribo', help='Path to ribo ORF file')
@click.option('--frame', help='Path to frame file')
@click.option('--prefix', help='Prefix to output file')
def theta_dist_cmd(rna, ribo, frame, prefix):
    theta_dist(rna, ribo, frame, prefix)


###################### theta rna distribution function #########################################
@cli.command(
    'theta-rna',
    context_settings=CONTEXT_SETTINGS,
    help='Test distribution of angles in RNA')
@click.option('--rna', help='Path to rna ORF file')
@click.option('--prefix', help='Prefix to output file')
@click.option('--cutoff', help='Cutoff of reads', type=int, default=10)
def theta_rna_cmd(rna, prefix, cutoff):
    theta_rna(rna, prefix, cutoff)


###################### parse-annotation function #########################################
@cli.command(
    'parse-annotation',
    context_settings=CONTEXT_SETTINGS,
    help='Parse annotation file to extract candidate ORFs')
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
    read_lengths = {
        28: 100023434,
        27: 2690000,
        30: 1234908,
        31: 790000,
        32: 8000,
        26: 123490
    }
    for length in read_lengths:
        metagenes[length] = s
    plot_metagene(metagenes, read_lengths, prefix, 100)
