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
import six
import pandas as pd

from . import __version__

from .count import export_gene_coverages
from .count import export_metagene_coverage
from .count import export_read_counts
from .count import merge_gene_coverages
from .count import merge_read_counts
from .count import export_read_length
from .count import read_enrichment
from .count import bedgraph_to_bigwig
from .count import bam_to_bedgraph
from .count import count_uniq_mapping_reads
from .count import extract_uniq_mapping_reads

from .orf import split_bam
from .orf import align_coverages
from .orf import merge_wigs
from .orf import prepare_orfs

from .infer_protocol import infer_protocol

from .sequence import export_gene_sequences

from .download import run_download_sra_script
from .coherence import naive_periodicity
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


###################### export-gene-coverages #################################
@cli.command(
    'export-gene-coverages',
    context_settings=CONTEXT_SETTINGS,
    help='Export gene level coverage for all genes of given region')
@click.option(
    '--bed',
    help='Path to bed file or a genome name (hg38_utr5, hg38_cds, hg38_utr3)',
    required=True)
@click.option('--pos_bw', help='Path to positive bigwig file', required=True)
@click.option('--neg_bw', help='Path to negative bigwig file', required=True)
@click.option(
    '--saveto', help='Path to write output', default=None, show_default=True)
@click.option(
    '--offset_5p',
    help='Number of upstream bases to count(5\')',
    type=int,
    default=0,
    show_default=True)
@click.option(
    '--offset_3p',
    help='Number of downstream bases to count(3\')',
    type=int,
    default=0,
    show_default=True)
def export_gene_coverages_cmd(bed, pos_bw, neg_bw, 
                              saveto, offset_5p, offset_3p):
    export_gene_coverages(bed, pos_bw, neg_bw, saveto, offset_5p, offset_3p)


###################### export-metagene-coverages ##############################
@cli.command(
    'export-metagene-coverage',
    context_settings=CONTEXT_SETTINGS,
    help='Export metagene coverage for given region')
@click.option(
    '--bed',
    help='Path to bed file or a genome name (hg38_utr5, hg38_cds, hg38_utr3)',
    required=True)
@click.option('--pos_bw', help='Path to positive bigwig file', required=True)
@click.option('--neg_bw', help='Path to negative bigwig file', required=True)
@click.option(
    '--max_positions',
    help='maximum positions to count',
    type=int,
    default=500,
    show_default=True)
@click.option(
    '--saveto', help='Path to write output', default=None, show_default=True)
@click.option(
    '--offset_5p',
    help='Number of upstream bases to count(5\')',
    type=int,
    default=0,
    show_default=True)
@click.option(
    '--offset_3p',
    help='Number of downstream bases to count(3\')',
    type=int,
    default=0,
    show_default=True)
def export_metagene_coverage_cmd(bed, pos_bw, neg_bw, max_positions, saveto,
                                 offset_5p, offset_3p):
    metagene_profile = export_metagene_coverage(bed, pos_bw, neg_bw,
                                                max_positions, saveto,
                                                offset_5p, offset_3p)

    for i, count in six.iteritems(metagene_profile):
        sys.stdout.write('{}\t{}'.format(i, count))
        sys.stdout.write(os.linesep)


###################### export-read-counts ##############################
@cli.command(
    'export-read-counts',
    context_settings=CONTEXT_SETTINGS,
    help='Export read counts from gene coverages file')
@click.option(
    '--gene_coverages', help='Path to gene coverages file', required=True)
@click.option(
    '--saveto', help='Path to write output', default=None, show_default=True)
@click.option(
    '--keep_offsets',
    help='whether keep the 5\' and 3\' offsets',
    is_flag=True)
def export_read_counts_cmd(gene_coverages, saveto, keep_offsets):
    export_read_counts(gene_coverages, saveto, keep_offsets)


###################### merge-gene-coverages ##############################
@cli.command(
    'merge-gene-coverages',
    context_settings=CONTEXT_SETTINGS,
    help='merge gene coverages to generate metagene coverage')
@click.option(
    '--gene_coverages', help='Path to gene coverages file', required=True)
@click.option(
    '--max_positions',
    help='maximum positions to count',
    type=int,
    default=500,
    show_default=True)
@click.option(
    '--saveto', help='Path to write output', default=None, show_default=True)
def merge_gene_coverages_cmd(gene_coverages, max_positions, saveto):
    merge_gene_coverages(gene_coverages, max_positions, saveto)


###################### merge-read-counts ##############################
@cli.command(
    'merge-read-counts',
    context_settings=CONTEXT_SETTINGS,
    help='merge read counts to generate count table')
@click.option(
    '--read_counts',
    help='Path to file containing read counts paths',
    required=True)
@click.option('--saveto', help='Path to write output', required=True)
def merge_read_counts_cmd(read_counts, saveto):
    merge_read_counts(read_counts, saveto)


#################### export-read-length ######################################
@cli.command(
    'export-read-length',
    context_settings=CONTEXT_SETTINGS,
    help='Calculate read length distribution')
@click.option('--bam', help='Path to BAM file', required=True)
@click.option('--saveto', help='Path to write read length dist tsv output')
def export_read_length_cmd(bam, saveto):
    counts = export_read_length(bam, saveto)
    for i, count in six.iteritems(dict(counts)):
        sys.stdout.write('{}\t{}'.format(i, count))
        sys.stdout.write(os.linesep)


#################### read-enrichment ######################################
@cli.command(
    'read-enrichment',
    context_settings=CONTEXT_SETTINGS,
    help='Calculate read enrichment for a certain range of lengths')
@click.option('--read_lengths', help='Path to read length tsv', required=True)
@click.option('--min_length', help='The low end of the range', default=28)
@click.option('--max_length', help='The high end of the range', default=32)
def read_enrichment_cmd(read_lengths, min_length, max_length):
    ratio = read_enrichment(read_lengths, min_length, max_length)
    sys.stdout.write('Enrichment of length range {}-{}: {}'.format(
        min_length, max_length, ratio))
    sys.stdout.write(os.linesep)


#################### periodicity #############################################
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


#################### export-gene-sequences ####################################
@cli.command(
    'export-gene-sequences',
    context_settings=CONTEXT_SETTINGS,
    help='Export gene level sequence for all genes of given region')
@click.option(
    '--bed',
    help='Path to bed file or a genome name (hg38_utr5, hg38_cds, hg38_utr3)',
    required=True)
@click.option('--fasta', help='Path to fasta file', required=True)
@click.option(
    '--saveto', help='Path to write output', default=None, show_default=True)
@click.option(
    '--offset_5p',
    help='Number of upstream bases to count(5\')',
    type=int,
    default=0)
@click.option(
    '--offset_3p',
    help='Number of downstream bases to count(3\')',
    type=int,
    default=0)
def export_gene_sequences_cmd(bed, fasta, saveto, offset_5p, offset_3p):
    export_gene_sequences(bed, fasta, saveto, offset_5p, offset_3p)


#################### split bam ####################################
@cli.command(
    'split-bam',
    context_settings=CONTEXT_SETTINGS,
    help='Split bam file by length and strand into wig files')
@click.option(
    '--bam',
    help='Path to bam file',
    required=True)
@click.option(
    '--protocol',
    help='Experimenta protocol [forward, reverse]',
    required=True)
@click.option(
    '--prefix', help='Prefix for output wig files')
def split_bam_cmd(bam, protocol, prefix):
    split_bam(bam, protocol, prefix)


#################### align coverages ####################################
@cli.command(
    'align-coverages',
    context_settings=CONTEXT_SETTINGS,
    help='Align coverages to a base to get relative lags')
@click.option(
    '--coverages',
    help='Path to file containing path to coverages',
    required=True)
@click.option(
    '--base',
    help='The reference length, usually the most abundant one',
    type=int,
    required=True)
@click.option(
    '--saveto', help='path to output tsv file')
def align_coverages_cmd(coverages, base, saveto):
    align_coverages(coverages, base, saveto)


#################### prepare orfs ####################################
@cli.command(
    'prepare-orfs',
    context_settings=CONTEXT_SETTINGS,
    help='Create bed file for all putative ORFs')
@click.option(
    '--gtf',
    help='Path to annotation file',
    required=True)
@click.option(
    '--fasta',
    help='Path to reference genome',
    required=True)
@click.option(
    '--prefix', help='Prefix to all output files')
def prepare_orfs_cmd(gtf, fasta, prefix):
    prepare_orfs(gtf, fasta, prefix)

#################### merge wigs ####################################
@cli.command(
    'merge-wigs',
    context_settings=CONTEXT_SETTINGS,
    help='merge wigs from all lengths by shifting offsets')
@click.option(
    '--wigs',
    help='Path to file containing path to wigs',
    required=True)
@click.option(
    '--offsets',
    help='Path to file containing offsets',
    required=True)
@click.option(
    '--strand',
    help='Strand of the wig files',
    required=True)
@click.option(
    '--saveto', help='path to merged wig file')
def merge_wigs_cmd(wigs, offsets, strand, saveto):
    merge_wigs(wigs, offsets, strand, saveto)


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
    default=None)
@click.option('--color', help='Color', default='royalblue')
@click.option(
    '--millify_labels',
    help='Convert labels on Y-axis to concise form?',
    is_flag=True)
@click.option('--identify_peak', help='Identify Peak?', is_flag=True)
@click.option(
    '--positions', help='Range of positions to plot', default='-60:500')
@click.option(
    '--saveto',
    help='Path to file (png/pdf) to save to',
    default=None,
    show_default=True)
@click.option(
    '--ylabel', help='Y axix label', default='Normalized RPF density')
@click.option(
    '--yrotation', default=0, help='Rotate y axis labels by', type=int)
@click.option(
    '--xrotation', default=0, help='Rotate x axis labels by', type=int)
@click.option('--ascii', help='Plot ascii', is_flag=True)
def plot_read_counts_cmd(counts, title, marker, color, millify_labels,
                         identify_peak, positions, saveto, ylabel, xrotation,
                         yrotation, ascii):
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
            ascii=ascii,
            xrotation=xrotation,
            yrotation=yrotation)
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
            input_is_stream=True,
            xrotation=xrotation,
            yrotation=yrotation)


######################## plot-read-length ####################################
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
    count = count_uniq_mapping_reads(bam)
    sys.stdout.write(str(count))
    sys.stdout.write(os.linesep)


###################### download function #########################################
@cli.command(
    'download-sra',
    context_settings=CONTEXT_SETTINGS,
    help='Download SRA data')
@click.option('--out', help='root directory to download all datasets')
@click.option('-a', '--ascp', help='Path to ascp private key')
@click.option(
    '-f',
    '--file',
    '--srpfile',
    help='File containing list of SRPs one per line')
@click.argument('srp_id_list', nargs=-1, required=True)
def download_srp_cmd(out, ascp, srpfile, srp_id_list):
    run_download_sra_script(out, ascp, srpfile, list(srp_id_list))


###################### infer-protocol function #########################################
@cli.command(
    'infer-protocol',
    context_settings=CONTEXT_SETTINGS,
    help='Infer protocol from BAM')
@click.option('--bam', help='Path to bam file')
@click.option(
    '--refseq',
    help='Path to reseq file to be used for defining the gene strands')
@click.option(
    '--prefix',
    help='Prefix to output protocol file')
@click.option(
    '--n_reads',
    type=int,
    default=20000,
    help='Number of mapped reads to use for estimation')
def infer_protocol_cmd(bam, refseq, prefix, n_reads):
    protocol, forward_mapped, reverse_mapped = infer_protocol(
        bam, refseq, prefix, n_reads)
    print(
        dedent('''\
                 Forward mapped proportion: {:.4f}
                 Reverse mapped proportion: {:.4f}
                 Likely protocol: {}'''.format(forward_mapped, reverse_mapped,
                                               protocol)))
