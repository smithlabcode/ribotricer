from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os
import sys

import click
from click_help_colors import HelpColorsGroup
import six
import pandas as pd

from . import __version__

from .coherence import naive_periodicity

from .count import bam_to_bedgraph
from .count import bedgraph_to_bigwig
from .count import unique_mapping_reads_count
from .count import extract_uniq_mapping_reads
from .count import export_gene_coverages
from .count import export_metagene_coverage
from .count import read_length_distribution

from .fasta import export_all_fasta

from .plotting import plot_read_counts
from .plotting import plot_read_length_dist

click.disable_unicode_literals_warning = True
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(
    cls=HelpColorsGroup,
    help_headers_color='yellow',
    help_options_color='green')
@click.version_option(version=__version__)
def cli():
    """riboraptor: Tool for ribosome profiling analysis"""
    pass


##############################################################################
#################### COUNT RELATED FUNCTIONS #################################
##############################################################################


###################### export-gene-coverages #################################
@cli.command(
    'export-gene-coverages',
    context_settings=CONTEXT_SETTINGS,
    help='Export gene level coverage for all genes for given region')
@click.option('--bigwig', '-bw', help='Path to bigwig', required=True)
@click.option('--region_bed', help='Path to bed file', required=True)
@click.option('--saveto', help='Path to write gene coverages tsv file')
@click.option(
    '--offset_5p',
    help='Number of upstream bases to count(5\')',
    type=int,
    default=60)
@click.option(
    '--offset_3p',
    help='Number of downstream bases to count(3\')',
    type=int,
    default=0)
@click.option(
    '--ignore_tx_version',
    help='Ignore version (.xyz) in gene names',
    is_flag=True)
def export_gene_coverages_cmd(bigwig, region_bed, saveto, offset_5p, offset_3p,
                              ignore_tx_version):
    export_gene_coverages(bigwig, region_bed, saveto, offset_5p, offset_3p,
                          ignore_tx_version)


###################### metagene-coverages #################################
@cli.command(
    'export-metagene-coverage',
    context_settings=CONTEXT_SETTINGS,
    help='Export metagene coverage for given region')
@click.option('--bigwig', '-bw', help='Path to bigwig', required=True)
@click.option(
    '--region_bed',
    help='Path to bed file or a genome name (hg38_utr5, hg38_cds, hg38_utr3)',
    required=True)
@click.option('--saveto', help='Path to write metagene coverage tsv file')
@click.option(
    '--max_positions', help='maximum positions to count', default=500)
@click.option(
    '--offset_5p',
    help='Number of upstream bases to count(5\')',
    type=int,
    default=60)
@click.option(
    '--offset_3p',
    help='Number of downstream bases to count(3\')',
    type=int,
    default=0)
@click.option(
    '--ignore_tx_version',
    help='Ignore version (.xyz) in gene names',
    is_flag=True)
def export_metagene_coverage_cmd(bigwig, region_bed, max_positions, saveto,
                                 offset_5p, offset_3p, ignore_tx_version):
    metagene_profile = export_metagene_coverage(
        bigwig=bigwig,
        region_bed_f=region_bed,
        max_positions=max_positions,
        saveto=saveto,
        offset_5p=offset_5p,
        offset_3p=offset_3p,
        ignore_tx_version=ignore_tx_version)
    for i, count in six.iteritems(metagene_profile):
        sys.stdout.write('{}\t{}'.format(i, count))
        sys.stdout.write(os.linesep)


#################### periodicity #####################################
@cli.command(
    'periodicity',
    context_settings=CONTEXT_SETTINGS,
    help='Calculate periodicity')
@click.option('--counts', help='Path to counts file (if not stdin)')
def periodicity_cmd(counts):
    if counts:
        counts = pd.read_table(counts)
        counts = pd.Series(
            counts['count'].tolist(), index=counts['position'].tolist())
        periodicity = naive_periodicity(counts)
    else:
        periodicity = naive_periodicity(
            sys.stdin.readlines(), input_is_stream=True)
    sys.stdout.write('Periodicity: {}'.format(periodicity))
    sys.stdout.write(os.linesep)


@cli.command(
    'read-length-dist',
    context_settings=CONTEXT_SETTINGS,
    help='Calculate read length distribution')
@click.option('--bam', help='Path to BAM file', required=True)
@click.option('--saveto', help='Path to write read length dist tsv output')
def rld_cmd(bam, saveto):
    counts = read_length_distribution(bam, saveto)
    for i, count in six.iteritems(dict(counts)):
        sys.stdout.write('{}\t{}'.format(i, count))
        sys.stdout.write(os.linesep)


@cli.command(
    'export-bed-fasta',
    context_settings=CONTEXT_SETTINGS,
    help='Export gene level fasta from specified bed regions')
@click.option('--region_bed', help='Path to bed file', required=True)
@click.option('--fasta', help='Path to fasta file', required=True)
@click.option(
    '--prefix',
    '-o',
    help='Path to write output',
    default=None,
    show_default=True)
@click.option('--chrom_sizes', help='Path to chrom.sizes', required=True)
@click.option(
    '--offset_5p',
    help='Number of upstream bases to count(5\')',
    type=int,
    default=60)
@click.option(
    '--offset_3p',
    help='Number of downstream bases to count(3\')',
    type=int,
    default=0)
@click.option(
    '--ignore_tx_version',
    help='Ignore version (.xyz) in gene names',
    is_flag=True)
def export_all_fasta_cmd(region_bed, chrom_sizes, fasta, prefix, offset_5p,
                         offset_3p, ignore_tx_version):
    export_all_fasta(region_bed, chrom_sizes, fasta, prefix, offset_5p,
                     offset_3p, ignore_tx_version)


##############################################################################
#################### PLOTING RELATED FUNCTIONS ###############################
##############################################################################


######################## plot-metagene ####################################
@cli.command(
    'plot-metagene',
    context_settings=CONTEXT_SETTINGS,
    help='Plot metagene profile')
@click.option('--counts', help='Path to counts file (if not stdin)')
@click.option('--title', help='Plot Title')
@click.option(
    '--marker',
    help='Marker (o/x) for plots',
    type=click.Choice(['o', 'x']),
    default='o')
@click.option('--color', help='Color', default='royalblue')
@click.option(
    '--millify_labels',
    help='Convert labels on Y-axis to concise form?',
    is_flag=True)
@click.option('--identify_peak', help='Identify Peak?', is_flag=True)
@click.option(
    '--positions', help='Range of positions to plot', default='-60:100')
@click.option(
    '--saveto',
    help='Path to file (png/pdf) to save to',
    default=None,
    show_default=True)
@click.option(
    '--ylabel', help='Y axix label', default='Normalized RPF density')
@click.option('--ascii', help='Do not plot ascii', is_flag=True)
def plot_read_counts_cmd(counts, title, marker, color, millify_labels,
                         identify_peak, positions, saveto, ylabel, ascii):
    if counts:
        plot_read_counts(
            counts,
            title=title,
            marker=marker,
            color=color,
            millify_labels=millify_labels,
            position_range=positions,
            identify_peak=identify_peak,
            ylabel=ylabel,
            saveto=saveto,
            ascii=ascii)
    else:
        plot_read_counts(
            sys.stdin.readlines(),
            title=title,
            marker=marker,
            color=color,
            millify_labels=millify_labels,
            identify_peak=identify_peak,
            position_range=positions,
            saveto=saveto,
            ylabel=ylabel,
            ascii=ascii,
            input_is_stream=True)


######################## plot-read-dist ####################################
@cli.command(
    'plot-read-length',
    context_settings=CONTEXT_SETTINGS,
    help='Plot read length distribution')
@click.option('--read-lengths', help='Path to read length pickle file')
@click.option('--title', help='Plot Title')
@click.option(
    '--millify_labels',
    help='Convert labels on Y-axis to concise form?',
    is_flag=True)
@click.option(
    '--saveto',
    help='Path to file (png/pdf) to save to',
    default=None,
    show_default=True)
@click.option('--ascii', help='Do not plot ascii', is_flag=True)
def plot_read_length_dist_cmd(read_lengths, title, millify_labels, saveto,
                              ascii):
    if read_lengths:
        plot_read_length_dist(
            read_lengths,
            title=title,
            millify_labels=millify_labels,
            input_is_stream=False,
            saveto=saveto,
            ascii=ascii)
    else:
        plot_read_length_dist(
            sys.stdin.readlines(),
            title=title,
            millify_labels=millify_labels,
            input_is_stream=True,
            saveto=saveto,
            ascii=ascii)


##############################################################################
#################### BAMTOOLS STYLE FUNCTIONS ################################
##############################################################################


####################### bam-to-bedgraph ######################################
@cli.command(
    'bam-to-bedgraph',
    context_settings=CONTEXT_SETTINGS,
    help='Convert bam to bedgraph')
@click.option('--bam', '-i', help='Path to BAM file', required=True)
@click.option(
    '--strand',
    '-s',
    help='Count from strand of this type only',
    type=click.Choice(['+', '-', 'both']),
    default='both')
@click.option(
    '--end_type',
    '-e',
    help='Pileup 5\' / 3\'/ either ends',
    type=click.Choice(['5prime', '3prime', 'either']),
    default='5prime')
@click.option(
    '--saveto',
    '-o',
    help='Path to write bedgraph output',
    default=None,
    show_default=True)
def bam_to_bedgraph_cmd(bam, strand, end_type, saveto):
    bedgraph = bam_to_bedgraph(bam, strand, end_type, saveto)
    if saveto is None:
        sys.stdout.write(str(bedgraph))
        sys.stdout.write(os.linesep)


####################### uniq-bam ##############################################
@cli.command(
    'uniq-bam',
    context_settings=CONTEXT_SETTINGS,
    help='Create a new bam with unique mapping reads only')
@click.option('--inbam', required=True)
@click.option('--outbam', required=True)
def extract_uniq_mapping_reads_cmd(inbam, outbam):
    extract_uniq_mapping_reads(inbam, outbam)


####################### bedgraph-to-bigwig ######################################
@cli.command(
    'bedgraph-to-bigwig',
    context_settings=CONTEXT_SETTINGS,
    help='Convert bedgraph to bigwig')
@click.option(
    '--bedgraph',
    '-bg',
    '-i',
    help='Path to bedgraph file (optional)',
    default=None,
    show_default=True)
@click.option(
    '--sizes', '-s', help='Path to genome chrom.sizes file', required=True)
@click.option(
    '--saveto', '-o', help='Path to write bigwig output', required=True)
def bedgraph_to_bigwig_cmd(bedgraph, sizes, saveto):
    if bedgraph:
        bedgraph_to_bigwig(bedgraph, sizes, saveto)
    else:
        bedgraph_to_bigwig(sys.stdin.readlines(), sizes, saveto, True)


###################### uniq-mapping-count ######################################
@cli.command(
    'uniq-mapping-count',
    context_settings=CONTEXT_SETTINGS,
    help='Count number of unique mapping reads')
@click.option('--bam', help='Path to BAM file', required=True)
def uniq_mapping_cmd(bam):
    count = unique_mapping_reads_count(bam)
    sys.stdout.write(str(count))
    sys.stdout.write(os.linesep)
