from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import os
import sys
import click
import six

from . import bam_to_bedgraph
from . import bedgraph_to_bigwig
from . import read_enrichment
from . import gene_coverage
# from .utils import htseq_to_cpm
from . import mapping_reads_summary
from . import read_length_distribution

from .plotting import plot_read_counts
from .plotting import plot_read_length_dist
from click_help_colors import HelpColorsGroup

click.disable_unicode_literals_warning = True
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(cls=HelpColorsGroup,
             help_headers_color='yellow',
             help_options_color='green')
@click.version_option(version='0.1.0')
def cli():
    """riboraptor: Tool for ribosome profiling analysis"""
    pass


@cli.command('bam-to-bedgraph', context_settings=CONTEXT_SETTINGS,
             help='Convert bam to bedgraph')
@click.option('--bam',
              help='Path to BAM file',
              required=True)
@click.option('--strand',
              help='Count from strand of this type only',
              type=click.Choice(['+', '-', 'both']),
              default='both')
@click.option('--end_type',
              help='Pileup 5\' / 3\'/ either ends',
              type=click.Choice(['5prime', '3prime', 'either']),
              default='5prime')
@click.option('--saveto',
              help='Path to write bedgraph output',
              default=None)
def bam_to_bedgraph_cmd(bam, strand,
                        end_type, saveto):
    bedgraph = bam_to_bedgraph(bam, strand,
                               end_type, saveto)
    if saveto is None:
        sys.stdout.write(bedgraph)
        sys.stdout.write(os.linesep)


@cli.command('bedgraph-to-bigwig', context_settings=CONTEXT_SETTINGS,
             help='Convert bedgraph to bigwig')
@click.option('--bedgraph', '--bg',
              help='Path to bedgraph file (optional)',
              default=None)
@click.option('--chrom_sizes',
              help='Path to genome chrom.sizes file',
              required=True)
@click.option('--saveto',
              help='Path to write bigwig output',
              required=True)
def bedgraph_to_bigwig_cmd(bedgraph, chrom_sizes, saveto):
    if bedgraph:
        bedgraph_to_bigwig(bedgraph, chrom_sizes, saveto)
    else:
        bedgraph_to_bigwig(sys.stdout.readlines(), chrom_sizes, saveto, True)


@cli.command('gene-coverage', context_settings=CONTEXT_SETTINGS,
             help='Calculate coverage across a gene')
@click.option('--gene',
              help='Gene name',
              required=True)
@click.option('--bed',
              help='BED file with \'gene\' annotated in name column',
              required=True)
@click.option('--bigwig', '--bw',
              help='Path to bigwig',
              required=True)
@click.option('--offset',
              help='Number of upstream bases to count',
              type=int,
              default=0)
def gene_coverage_cmd(gene, bed, bigwig, offset):
    coverage_combined, _, _, gene_offset = gene_coverage(gene, bed,
                                                         bigwig, offset)
    for l, count in six.iteritems(dict(coverage_combined)):
        sys.stdout.write('{}\t{}'.format(l, count))
        sys.stdout.write(os.linesep)


@cli.command('mapping-summary', context_settings=CONTEXT_SETTINGS,
             help='Mapping summary')
@click.option('--bam',
              help='Path to BAM file',
              required=True)
def mapping_reads_summary_cmd(bam):
    counts = mapping_reads_summary(bam)
    for l, count in six.iteritems(dict(counts)):
        sys.stdout.write('{}\t{}'.format(l, count))
        sys.stdout.write(os.linesep)


@cli.command('plot-read-counts', context_settings=CONTEXT_SETTINGS,
             help='Plot read counts distribution across a gene')
@click.option('--millify_labels',
              help='Convert labels on Y-axis to concise form?',
              is_flag=True)
@click.option('--identify_peak',
              help='Identify Peak?',
              is_flag=True)
@click.option('--marker',
              help='Marker (o/x) for plots',
              type=click.Choice(['o', 'x']),
              default='o')
@click.option('--color',
              help='Color',
              default='royalblue')
@click.option('--saveto',
              help='Path to file (png/pdf) to save to',
              required=True)
def plot_read_counts_cmd(counts,
                         marker, color,
                         millify_labels,
                         identify_peak,
                         saveto):
    plot_read_counts(counts,
                     marker=marker, color=color,
                     millify_labels=millify_labels,
                     identify_peak=identify_peak, saveto=saveto)


@cli.command('plot-read-dist', context_settings=CONTEXT_SETTINGS,
             help='Plot read length distribution')
@click.option('--millify_labels',
              help='Convert labels on Y-axis to concise form?',
              is_flag=True)
@click.option('--saveto',
              help='Path to file (png/pdf) to save to',
              required=True)
def plot_read_length_dist_cmd(millify_labels, saveto):
    plot_read_length_dist(sys.stdin.readlines(), millify_labels=millify_labels,
                          input_is_stream=True, saveto=saveto)


@cli.command('read-enrichment', context_settings=CONTEXT_SETTINGS,
             help='Calculate read length enrichment')
@click.option('--lrange',
              help='reads lengths to use for enrichment',
              default='28-32')
def read_enrichment_cmd(lrange):
    enrichment, pvalue = read_enrichment(sys.stdin.readlines(),
                                         lrange,
                                         True)
    sys.stdout.write('(Enrichment: {}, pval: {})'.format(enrichment, pvalue))
    sys.stdout.write(os.linesep)


@cli.command('read-length-dist', context_settings=CONTEXT_SETTINGS,
             help='Calculate read length distribution')
@click.option('--bam',
              help='Path to BAM file',
              required=True)
def rld_cmd(bam):
    counts = read_length_distribution(bam)
    for l, count in six.iteritems(dict(counts)):
        sys.stdout.write('{}\t{}'.format(l, count))
        sys.stdout.write(os.linesep)
