from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os
import sys

import click
from click_help_colors import HelpColorsGroup
import six
import pandas as pd

from . import __version__
from .coherence import get_periodicity
from .count import bam_to_bedgraph
from .count import bedgraph_to_bigwig
from .count import count_reads_per_gene
from .count import count_reads_per_gene_htseq
from .count import collapse_gene_coverage_to_metagene
from .count import count_utr5_utr3_cds
from .count import diff_region_enrichment
from .count import export_gene_coverages
from .count import export_single_gene_coverage
from .count import extract_uniq_mapping_reads
from .count import gene_coverage
from .count import htseq_to_tpm
from .count import mapping_reads_summary
from .count import metagene_coverage
from .count import pickle_bed_file
from .count import read_enrichment
from .count import read_length_distribution
from .count import unique_mapping_reads_count

from .fasta import export_all_fasta
from .fasta import complete_gene_fasta

from .helpers import parse_star_logs

from .plotting import plot_framewise_counts
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


@cli.command(
    'uniq-bam',
    context_settings=CONTEXT_SETTINGS,
    help='Create a new bam with unique mapping reads only')
@click.option('--inbam', required=True)
@click.option('--outbam', required=True)
def extract_uniq_mapping_reads_cmd(inbam, outbam):
    extract_uniq_mapping_reads(inbam, outbam)


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


@cli.command(
    'export-gene-coverages',
    context_settings=CONTEXT_SETTINGS,
    help='Export gene coverages')
@click.option('--bigwig', '-bw', help='Path to bigwig', required=True)
@click.option('--region_bed', help='Path to CDS bed file', required=True)
@click.option('--prefix', help='Save gene coverages file to')
@click.option(
    '--offset_5p',
    help='Number of upstream bases to count(5\')',
    type=int,
    default=60)
@click.option(
    '--offset_3p',
    help='Number of upstream bases to count(3\')',
    type=int,
    default=0)
@click.option(
    '--ignore_tx_version',
    help='Ignore version (.xyz) in gene names',
    is_flag=True)
def export_gene_coverages_cmd(bigwig, region_bed, prefix, offset_5p, offset_3p,
                              ignore_tx_version):
    export_gene_coverages(bigwig, region_bed, prefix, offset_5p, offset_3p,
                          ignore_tx_version)


@cli.command(
    'export-single-gene-coverage',
    context_settings=CONTEXT_SETTINGS,
    help='Export coverage for a gene')
@click.option('--bigwig', '-bw', help='Path to bigwig', required=True)
@click.option('--region_bed', help='Path to CDS bed file', required=True)
@click.option('--gene', help='Gene id', required=True)
@click.option('--prefix', help='Save gene coverages file to', required=True)
@click.option(
    '--offset_5p',
    help='Number of upstream bases to count(5\')',
    type=int,
    default=60)
@click.option(
    '--offset_3p',
    help='Number of upstream bases to count(3\')',
    type=int,
    default=0)
@click.option(
    '--ignore_tx_version',
    help='Ignore version (.xyz) in gene names',
    is_flag=True)
def export_single_gene_coverage_cmd(bigwig, region_bed, gene, prefix,
                                    offset_5p, offset_3p, ignore_tx_version):
    export_single_gene_coverage(bigwig, region_bed, gene, prefix, offset_5p,
                                offset_3p, ignore_tx_version)


@cli.command(
    'mapping-summary',
    context_settings=CONTEXT_SETTINGS,
    help='Mapping summary')
@click.option('--bam', help='Path to BAM file', required=True)
@click.option('--prefix', help='Save pickle files to')
def mapping_reads_summary_cmd(bam, prefix):
    counts = mapping_reads_summary(bam, prefix)
    for l, count in six.iteritems(dict(counts)):
        sys.stdout.write('{}\t{}'.format(l, count))
        sys.stdout.write(os.linesep)


@cli.command(
    'metagene-coverage',
    context_settings=CONTEXT_SETTINGS,
    help='Calculate metagene coverage')
@click.option('--bigwig', '-bw', help='Path to bigwig', required=True)
@click.option('--region_bed', help='Path to CDS bed file', required=True)
@click.option(
    '--max-positions',
    help='Number of upstream bases to count',
    type=int,
    default=200)
@click.option('--htseq_f', help='Path to htseq counts file')
@click.option('--prefix', help='Save pickle files to')
@click.option(
    '--offset_5p',
    help='Number of upstream bases to count(5\')',
    type=int,
    default=60)
@click.option(
    '--offset_3p',
    help='Number of upstream bases to count(3\')',
    type=int,
    default=0)
@click.option(
    '--n-meta',
    help='Number of genes to use for calculating metagene average',
    type=int,
    default=-1)
@click.option(
    '--n-save-gene',
    help='Number of genes profiles to pickle',
    type=int,
    default=0)
@click.option(
    '--ignore_tx_version',
    help='Ignore version (.xyz) in gene names',
    is_flag=True)
def metagene_coverage_cmd(bigwig, region_bed, max_positions, htseq_f, prefix,
                          offset_5p, offset_3p, n_meta, n_save_gene,
                          ignore_tx_version):
    metagene_profile = metagene_coverage(
        bigwig=bigwig,
        region_bed_f=region_bed,
        max_positions=max_positions,
        htseq_f=htseq_f,
        prefix=prefix,
        offset_5p=offset_5p,
        offset_3p=offset_3p,
        top_n_meta=n_meta,
        top_n_gene=n_save_gene,
        ignore_tx_version=ignore_tx_version)
    for l, count in six.iteritems(metagene_profile):
        sys.stdout.write('{}\t{}'.format(l, count))
        sys.stdout.write(os.linesep)


@cli.command(
    'periodicity',
    context_settings=CONTEXT_SETTINGS,
    help='Calculate periodicity')
@click.option('--counts', help='Path to counts file (if not stdin)')
def periodicity_cmd(counts):
    if counts:
        periodicity, pval = get_periodicity(counts)
    else:
        periodicity, pval = get_periodicity(
            sys.stdin.readlines(), input_is_stream=True)
    sys.stdout.write('Periodicity: {}({})'.format(periodicity, pval))
    sys.stdout.write(os.linesep)


@cli.command(
    'plot-read-counts',
    context_settings=CONTEXT_SETTINGS,
    help='Plot read counts distribution across a gene')
@click.option('--counts', help='Path to counts file (if not stdin)')
@click.option('--title', help='Plot Title', required=True)
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
@click.option('--ylabel', help='Y axix label')
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


@cli.command(
    'plot-framewise-counts',
    context_settings=CONTEXT_SETTINGS,
    help='Plot read counts highlighting frames')
@click.option('--counts', help='Path to counts file (if not stdin)')
@click.option('--title', help='Plot Title', required=True)
@click.option(
    '--millify_labels',
    help='Convert labels on Y-axis to concise form?',
    is_flag=True)
@click.option(
    '--positions', help='Range of positions to plot', default='-60:100')
@click.option(
    '--saveto',
    help='Path to file (png/pdf) to save to',
    default=None,
    show_default=True)
@click.option('--ascii', help='Plot ascii', is_flag=True)
def plot_framewise_counts_cmd(counts, title, millify_labels, positions, saveto,
                              ascii):
    if counts:
        plot_framewise_counts(
            counts,
            title=title,
            millify_labels=millify_labels,
            position_range=positions,
            saveto=saveto,
            ascii=ascii)
    else:
        plot_framewise_counts(
            sys.stdin.readlines(),
            title=title,
            millify_labels=millify_labels,
            position_range=positions,
            saveto=saveto,
            ascii=ascii,
            input_is_stream=True)


@cli.command(
    'plot-read-dist',
    context_settings=CONTEXT_SETTINGS,
    help='Plot read length distribution')
@click.option('--read-lengths', help='Path to read length pickle file')
@click.option('--title', help='Plot Title', required=True)
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


@cli.command(
    'read-enrichment',
    context_settings=CONTEXT_SETTINGS,
    help='Calculate read length enrichment')
@click.option(
    '--lrange', help='reads lengths to use for enrichment', default='28-32')
@click.option(
    '--infile', '-i', help='Tab separated read length distribution file')
def read_enrichment_cmd(infile, lrange):
    if not infile:
        enrichment, pvalue = read_enrichment(sys.stdin.readlines(), lrange,
                                             True)
    else:
        enrichment, pvalue = read_enrichment(infile, lrange, False, True)

    sys.stdout.write('(Enrichment: {}, pval: {})'.format(enrichment, pvalue))
    sys.stdout.write(os.linesep)


@cli.command(
    'read-length-dist',
    context_settings=CONTEXT_SETTINGS,
    help='Calculate read length distribution')
@click.option('--bam', help='Path to BAM file', required=True)
@click.option('--prefix', help='Prefix to write pickled contents')
def rld_cmd(bam, prefix):
    counts = read_length_distribution(bam, prefix)
    for l, count in six.iteritems(dict(counts)):
        sys.stdout.write('{}\t{}'.format(l, count))
        sys.stdout.write(os.linesep)


@cli.command(
    'uniq-mapping-count',
    context_settings=CONTEXT_SETTINGS,
    help='Count number of unique mapping reads')
@click.option('--bam', help='Path to BAM file', required=True)
def uniq_mapping_cmd(bam):
    count = unique_mapping_reads_count(bam)
    sys.stdout.write(str(count))
    sys.stdout.write(os.linesep)


@cli.command(
    'diff-region-enrichment',
    context_settings=CONTEXT_SETTINGS,
    help='calculate enrichment of cds over utr3/utr5 counts and alike')
@click.option(
    '--numerator', help='path to counts file for numerator', required=True)
@click.option(
    '--denominator', help='path to counts file for denominator', required=True)
@click.option('--prefix', help='prefix to write pickled contents')
def diff_region_enrichment_cmd(numerator, denominator, prefix):
    enrichment = diff_region_enrichment(numerator, denominator, prefix)
    for l, e in six.iteritems(dict(enrichment)):
        sys.stdout.write('{}\t{}'.format(l, e))
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


@cli.command(
    'export-complete-fasta',
    context_settings=CONTEXT_SETTINGS,
    help='Export gene level fasta from specified bed regions')
@click.option('--utr5_bed', help='Path to  UTR5 bed file', required=True)
@click.option('--cds_bed', help='Path to  UTR5 bed file', required=True)
@click.option('--utr3_bed', help='Path to  UTR5 bed file', required=True)
@click.option('--fasta', help='Path to fasta file', required=True)
@click.option(
    '--prefix',
    '-o',
    help='Path to write output',
    default=None,
    show_default=True)
def export_complete_fasta_cmd(utr5_bed, cds_bed, utr3_bed, fasta, prefix):
    complete_gene_fasta(utr5_bed, cds_bed, utr3_bed, fasta, prefix)


@cli.command(
    'extract-star-logs',
    context_settings=CONTEXT_SETTINGS,
    help='collpase star logs to a dataframe')
@click.option('--starlogs', help='Path to star.logs.out file', required=True)
@click.option('--outfile', help='Path to output dataframe.tsv', required=True)
def extract_star_logs(starlogs, outfile):
    parse_star_logs(starlogs, outfile)


@cli.command(
    'collapse-gene-coverage',
    context_settings=CONTEXT_SETTINGS,
    help='Collapse gene coverage to metagene of target length',
)
@click.option(
    '--gene_coverage', help='Path to gene coverage tsv', required=True)
@click.option(
    '--target_length',
    help='Target length of metagene',
    required=True,
    type=int)
@click.option('--outfile', help='Path to output file', required=True)
def collapse_gene_coverage_to_metagene_cmd(gene_coverage, target_length,
                                           outfile):
    collapse_gene_coverage_to_metagene(gene_coverage, target_length, outfile)
