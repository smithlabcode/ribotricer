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
from .count import count_reads_in_features
from .count import count_utr5_utr3_cds
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
@click.option('--bam', '-i',
              help='Path to BAM file',
              required=True)
@click.option('--strand', '-s',
              help='Count from strand of this type only',
              type=click.Choice(['+', '-', 'both']),
              default='both')
@click.option('--end_type', '-e',
              help='Pileup 5\' / 3\'/ either ends',
              type=click.Choice(['5prime', '3prime', 'either']),
              default='5prime')
@click.option('--saveto', '-o',
              help='Path to write bedgraph output',
              default=None, show_default=True)
def bam_to_bedgraph_cmd(bam, strand,
                        end_type, saveto):
    bedgraph = bam_to_bedgraph(bam, strand,
                               end_type, saveto)
    if saveto is None:
        sys.stdout.write(str(bedgraph))
        sys.stdout.write(os.linesep)


@cli.command('bedgraph-to-bigwig', context_settings=CONTEXT_SETTINGS,
             help='Convert bedgraph to bigwig')
@click.option('--bedgraph', '-bg', '-i',
              help='Path to bedgraph file (optional)',
              default=None, show_default=True)
@click.option('--sizes', '-s',
              help='Path to genome chrom.sizes file',
              required=True)
@click.option('--saveto', '-o',
              help='Path to write bigwig output',
              required=True)
def bedgraph_to_bigwig_cmd(bedgraph, sizes, saveto):
    if bedgraph:
        bedgraph_to_bigwig(bedgraph, sizes, saveto)
    else:
        bedgraph_to_bigwig(sys.stdin.readlines(), sizes, saveto, True)


@cli.command('count-in-feature', context_settings=CONTEXT_SETTINGS,
               help='Count reads in given feature bed file')
@click.option('--bam', help='Path to bam file', required=True)
@click.option('--bed', help='Path to 5\'utr file', required=True)
def count_reads_in_features_cmd(bam, bed):
    counts = count_reads_in_features(bam, bed)
    sys.stdout.write('{}'.format(counts))
    sys.stdout.write(os.linesep)


@cli.command('count-all-features', context_settings=CONTEXT_SETTINGS,
               help='Count reads in 5\'UTr/CDs/3\'UTR regions')
@click.option('--bam', help='Path to bam file', required=True)
@click.option('--utr5-bed', help='Path to 5\'utr file')
@click.option('--cds-bed', help='Path to CDS file')
@click.option('--utr3-bed', help='Path to 3\'UTR file')
@click.option('--genome', '-g', help='Genome (for loading bed files internally)')
@click.option('-s', help='Force strandedness', is_flag=True)
@click.option('--prefix', help='Prefix to write pickled contents')
def count_utr5_utr3_cds_cmd(bam, utr5_bed, cds_bed, utr3_bed,
                            genome, s, prefix):
    counts = count_utr5_utr3_cds(bam=bam, utr5_bed=utr5_bed,
                                 cds_bed=cds_bed, utr3_bed=utr3_bed,
                                 genome=genome, force_strandedness=s,
                                 saveto=prefix)
    for region, count in six.iteritems(dict(counts)):
        sys.stdout.write('{}\t{}'.format(region, count))
        sys.stdout.write(os.linesep)


@cli.command('gene-coverage', context_settings=CONTEXT_SETTINGS,
             help='Calculate coverage across a gene')
@click.option('--gene', '-n',
              help='Gene name',
              required=True)
@click.option('--bed', '-i',
              help='BED file with \'gene\' annotated in name column',
              required=True)
@click.option('--bigwig', '-bw',
              help='Path to bigwig',
              required=True)
@click.option('--offset', '-o',
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
@click.option('--bigwig', '-bw',
              help='Path to bigwig',
              required=True)
@click.option('--region_bed',
              help='Path to CDS bed file',
              required=True)
@click.option('--max-positions',
              help='Number of upstream bases to count',
              type=int,
              default=200)
@click.option('--htseq_f',
              help='Path to htseq counts file')
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
                          region_bed,
                          max_positions,
                          htseq_f,
                          prefix,
                          offset,
                          n_meta,
                          n_save_gene,
                          ignore_tx_version):
    metagene_profile = metagene_coverage(bigwig=bigwig,
                                         region_bed_f=region_bed,
                                         max_positions=max_positions,
                                         htseq_f=htseq_f,
                                         prefix=prefix,
                                         offset=offset,
                                         top_n_meta=n_meta,
                                         top_n_gene=n_save_gene,
                                         ignore_tx_version=ignore_tx_version)
    for l, count in metagene_profile.iteritems():
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
