from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from collections import OrderedDict
import os
import sys

import click
from click_help_colors import HelpColorsGroup
import six

from .count import bam_to_bedgraph
from .count import bedgraph_to_bigwig
from .count import read_enrichment
from .count import gene_coverage
# from .count import htseq_to_cpm
from .count import mapping_reads_summary
from .count import metagene_coverage
from .count import read_length_distribution
from .count import unique_mapping_reads_count

from .plotting import plot_read_counts
from .plotting import plot_read_length_dist

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


@cli.command('metagene-coverage', context_settings=CONTEXT_SETTINGS,
             help='Plot metagene plot')
@click.option('--bigwig', '--bw',
              help='Path to bigwig',
              required=True)
@click.option('--htseq_f',
              help='Path to htseq counts file',
              required=True)
@click.option('--region_bed',
              help='Path to CDS bed file',
              required=True)
@click.option('--prefix',
              help='Save pickle files to')
@click.option('--offset',
              help='Number of upstream bases to count',
              type=int,
              default=60)
@click.option('--n-meta',
              help='Number of genes to use for calculating metagene average',
              type=int,
              default=-1)
@click.option('--n-save-gene',
              help='Number of genes profiles to pickle',
              type=int,
              default=0)
@click.option('--ignore_tx_version',
              help='Ignore version (.xyz) in gene names',
              is_flag=True)
def metagene_coverage_cmd(bigwig,
                          htseq_f,
                          region_bed,
                          prefix,
                          offset,
                          n_meta,
                          n_save_gene,
                          ignore_tx_version):
    metagene_profile = metagene_coverage(bigwig,
                                         htseq_f,
                                         region_bed,
                                         prefix,
                                         master_offset=60,
                                         top_n_meta=-1,
                                         top_n_gene=10,
                                         ignore_tx_version=True)
    for l, count in six.iteritems(metagene_profile.to_dict(OrderedDict)):
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
@click.option('--infile', '-i', help='Tab separated read length distribution file')
def read_enrichment_cmd(infile, lrange):
    if not infile:
        enrichment, pvalue = read_enrichment(sys.stdin.readlines(),
                                             lrange,
                                             True)
    else:
        enrichment, pvalue = read_enrichment(infile,
                                             lrange,
                                             False,
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


@cli.command('uniq-mapping-count', context_settings=CONTEXT_SETTINGS,
             help='Count number of unique mapping reads')
@click.option('--bam',
              help='Path to BAM file',
              required=True)
def uniq_mapping_cmd(bam):
    count = unique_mapping_reads_count(bam)
    sys.stdout.write(str(count))
    sys.stdout.write(os.linesep)
