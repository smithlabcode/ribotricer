"""Command line interface for ribotricer
"""
# Part of ribotricer software
#
# Copyright (C) 2020 Saket Choudhary, Wenzheng Li, and Andrew D Smith
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

import click
import os
import sys

from . import __version__
from .common import _clean_input
from .const import CUTOFF
from .const import MINIMUM_VALID_CODONS
from .const import MINIMUM_VALID_CODONS_RATIO
from .const import MINIMUM_READS_PER_CODON
from .const import MINIMUM_DENSITY_OVER_ORF

from .count_orfs import count_orfs
from .count_orfs import count_orfs_codon
from .detect_orfs import detect_orfs
from .learn_cutoff import determine_cutoff_bam
from .learn_cutoff import determine_cutoff_tsv

from .orf_seq import orf_seq
from .prepare_orfs import prepare_orfs

from click_help_colors import HelpColorsGroup

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(
    cls=HelpColorsGroup, help_headers_color="yellow", help_options_color="green"
)
@click.version_option(version=__version__)
def cli():
    """ribotricer: Tool for detecting translating ORF from Ribo-seq data"""
    pass


###################### prepare-orfs function #########################################
@cli.command(
    "prepare-orfs",
    context_settings=CONTEXT_SETTINGS,
    help="Extract candidate ORFS based on GTF and FASTA files",
)
@click.option("--gtf", help="Path to GTF file", required=True)
@click.option("--fasta", help="Path to FASTA file", required=True)
@click.option("--prefix", help="Prefix to output file", required=True)
@click.option(
    "--min_orf_length",
    type=int,
    default=60,
    show_default=True,
    help="The minimum length (nts) of ORF to include",
)
@click.option(
    "--start_codons",
    default="ATG",
    show_default=True,
    help="Comma separated list of start codons",
)
@click.option(
    "--stop_codons",
    default="TAG,TAA,TGA",
    show_default=True,
    help="Comma separated list of stop codons",
)
@click.option(
    "--longest",
    help="Choose the most upstream start codon if multiple in frame ones exist",
    is_flag=True,
)
def prepare_orfs_cmd(
    gtf, fasta, prefix, min_orf_length, start_codons, stop_codons, longest
):
    if not os.path.isfile(gtf):
        sys.exit("Error: GTF file not found")

    if not os.path.isfile(fasta):
        sys.exit("Error: FASTA file not found")

    if min_orf_length <= 0:
        sys.exit("Error: min ORF length at least to be 1")

    start_codons = set([x.strip().upper() for x in start_codons.strip().split(",")])
    if not start_codons:
        sys.exit("Error: start codons cannot be empty")
    if not all([len(x) == 3 and set(x) <= {"A", "C", "G", "T"} for x in start_codons]):
        sys.exit("Error: invalid codon, only A, C, G, T allowed")

    stop_codons = set([x.strip().upper() for x in stop_codons.strip().split(",")])
    if not stop_codons:
        sys.exit("Error: stop codons cannot be empty")
    if not all([len(x) == 3 and set(x) <= {"A", "C", "G", "T"} for x in stop_codons]):
        sys.exit("Error: invalid codon, only A, C, G, T allowed")

    print("Using start codons: {}".format(",".join(start_codons)))
    prepare_orfs(gtf, fasta, prefix, min_orf_length, start_codons, stop_codons, longest)


###################### detect-orfs function #########################################
@cli.command(
    "detect-orfs",
    context_settings=CONTEXT_SETTINGS,
    help="Detect translating ORFs from BAM file",
)
@click.option("--bam", help="Path to BAM file", required=True)
@click.option(
    "--ribotricer_index",
    help=(
        "Path to the index file of ribotricer\n"
        "This file should be generated using ribotricer prepare-orfs"
    ),
    required=True,
)
@click.option("--prefix", help="Prefix to output file", required=True)
@click.option(
    "--stranded",
    type=click.Choice(["yes", "no", "reverse"]),
    default=None,
    show_default=True,
    help=(
        "whether the data is from a strand-specific assay"
        " If not provided, the experimental protocol will be automatically inferred"
    ),
)
@click.option(
    "--read_lengths",
    default=None,
    show_default=True,
    help=(
        "Comma separated read lengths to be used, such as 28,29,30\n"
        "If not provided, it will be automatically determined by assessing"
        " the metagene periodicity"
    ),
)
@click.option(
    "--psite_offsets",
    default=None,
    show_default=True,
    help=(
        "Comma separated P-site offsets for each read length "
        "matching the read lengths provided.\n"
        "If not provided, reads from different read lengths will be "
        "automatically aligned using cross-correlation"
    ),
)
@click.option(
    "--phase_score_cutoff",
    type=float,
    default=CUTOFF,
    show_default=True,
    help="Phase score cutoff for determining active translation",
)
@click.option(
    "--min_valid_codons",
    type=int,
    default=MINIMUM_VALID_CODONS,
    show_default=True,
    help="Minimum number of codons with non-zero reads for determining active translation",
)
@click.option(
    "--min_reads_per_codon",
    type=int,
    default=MINIMUM_READS_PER_CODON,
    show_default=True,
    help="Minimum number of reads per codon for determining active translation",
)
@click.option(
    "--min_valid_codons_ratio",
    type=float,
    default=MINIMUM_VALID_CODONS_RATIO,
    show_default=True,
    help="Minimum ratio of codons with non-zero reads to total codons for determining active translation",
)
@click.option(
    "--min_read_density",
    type=float,
    default=MINIMUM_DENSITY_OVER_ORF,
    show_default=True,
    help="Minimum read density (total_reads/length) over an ORF total codons for determining active translation",
)
@click.option(
    "--report_all",
    help=("Whether output all ORFs including those " "non-translating ones"),
    is_flag=True,
)
def detect_orfs_cmd(
    bam,
    ribotricer_index,
    prefix,
    stranded,
    read_lengths,
    psite_offsets,
    phase_score_cutoff,
    min_valid_codons,
    min_reads_per_codon,
    min_valid_codons_ratio,
    min_read_density,
    report_all,
):
    if not os.path.isfile(bam):
        sys.exit("Error: BAM file not found")

    if not os.path.isfile(ribotricer_index):
        sys.exit("Error: ribotricer index file not found")

    if read_lengths is not None:
        try:
            read_lengths = [int(x.strip()) for x in read_lengths.strip().split(",")]
        except:
            sys.exit("Error: cannot convert read_lengths into integers")
        if not all([x > 0 for x in read_lengths]):
            sys.exit("Error: read length must be positive")

    if read_lengths is None and psite_offsets is not None:
        sys.exit("Error: psite_offsets only allowed when read_lengths is provided")
    if read_lengths is not None and psite_offsets is not None:
        try:
            psite_offsets = [int(x.strip()) for x in psite_offsets.strip().split(",")]
        except:
            sys.exit("Error: cannot convert psite_offsets into integers")
        if len(read_lengths) != len(psite_offsets):
            sys.exit("Error: psite_offsets must match read_lengths")
        if not all(x >= 0 for x in psite_offsets):
            sys.exit("Error: P-site offset must be >= 0")
        if not all(x > y for (x, y) in zip(read_lengths, psite_offsets)):
            sys.exit("Error: P-site offset must be smaller than read length")
        psite_offsets = dict(list(zip(read_lengths, psite_offsets)))
    if stranded == "yes":
        stranded = "forward"
    detect_orfs(
        bam,
        ribotricer_index,
        prefix,
        stranded,
        read_lengths,
        psite_offsets,
        phase_score_cutoff,
        min_valid_codons,
        min_reads_per_codon,
        min_valid_codons_ratio,
        min_read_density,
        report_all,
    )


###################### count-orfs function #########################################
@cli.command(
    "count-orfs",
    context_settings=CONTEXT_SETTINGS,
    help="Count reads for detected ORFs at gene level",
)
@click.option(
    "--ribotricer_index",
    help=(
        "Path to the index file of ribotricer\n"
        "This file should be generated using ribotricer prepare-orfs"
    ),
    required=True,
)
@click.option(
    "--detected_orfs",
    help=(
        "Path to the detected orfs file\n"
        "This file should be generated using ribotricer detect-orfs"
    ),
    required=True,
)
@click.option("--features", help="ORF types separated with comma", required=True)
@click.option("--out", help="Path to output file", required=True)
@click.option(
    "--report_all",
    help=("Whether output all ORFs including those " "non-translating ones"),
    is_flag=True,
)
def count_orfs_cmd(ribotricer_index, detected_orfs, features, out, report_all):

    if not os.path.isfile(ribotricer_index):
        sys.exit("Error: ribotricer index file not found")

    if not os.path.isfile(detected_orfs):
        sys.exit("Error: detected orfs file not found")

    features = set(x.strip() for x in features.strip().split(","))

    count_orfs(ribotricer_index, detected_orfs, features, out, report_all)


###################### count-orfs-codon function #########################################
@cli.command(
    "count-orfs-codon",
    context_settings=CONTEXT_SETTINGS,
    help="Count reads for detected ORFs at codon level",
)
@click.option(
    "--ribotricer_index",
    help=(
        "Path to the index file of ribotricer\n"
        "This file should be generated using ribotricer prepare-orfs"
    ),
    required=True,
)
@click.option(
    "--detected_orfs",
    help=(
        "Path to the detected orfs file\n"
        "This file should be generated using ribotricer detect-orfs"
    ),
    required=True,
)
@click.option("--features", help="ORF types separated with comma", required=True)
@click.option("--ribotricer_index_fasta", help="Path to ORF seq file", required=True)
@click.option("--prefix", help="Prefix for output files", required=True)
@click.option(
    "--report_all",
    help=("Whether output all ORFs including those " "non-translating ones"),
    is_flag=True,
)
def count_orfs_codon_cmd(
    ribotricer_index,
    detected_orfs,
    features,
    ribotricer_index_fasta,
    prefix,
    report_all,
):

    if not os.path.isfile(ribotricer_index):
        sys.exit("Error: ribotricer index file not found")

    if not os.path.isfile(detected_orfs):
        sys.exit("Error: detected orfs file not found")

    if not os.path.isfile(ribotricer_index_fasta):
        sys.exit("Error: ribotricer_index_fasta file not found")

    features = set(x.strip() for x in features.strip().split(","))

    count_orfs_codon(
        ribotricer_index,
        detected_orfs,
        features,
        ribotricer_index_fasta,
        prefix,
        report_all,
    )


###################### orfs-seq function #########################################
@cli.command(
    "orfs-seq",
    context_settings=CONTEXT_SETTINGS,
    help="Generate sequence for ORFs in ribotricer's index",
)
@click.option(
    "--ribotricer_index",
    help=(
        "Path to the index file of ribotricer\n"
        "This file should be generated using ribotricer prepare-orfs"
    ),
    required=True,
)
@click.option("--fasta", help="Path to FASTA file", required=True)
@click.option(
    "--protein", help="Output protein sequence instead of nucleotide", is_flag=True
)
@click.option("--saveto", help="Path to output file", required=True)
def orf_seq_cmd(ribotricer_index, fasta, saveto, protein):
    if not os.path.isfile(ribotricer_index):
        sys.exit("Error: ribotricer index file not found")

    if not os.path.isfile(fasta):
        sys.exit("Error: fasta file not found")

    orf_seq(ribotricer_index, fasta, saveto, protein)


###################### learn-cutoff function #########################################
@cli.command(
    "learn-cutoff",
    context_settings=CONTEXT_SETTINGS,
    help="Learn phase score cutoff from BAM/TSV file",
)
@click.option("--ribo_bams", help="Path(s) to Ribo-seq BAM file separated by comma")
@click.option("--rna_bams", help="Path(s) to RNA-seq BAM file separated by comma")
@click.option(
    "--ribo_tsvs",
    help="Path(s) to Ribo-seq *_translating_ORFs.tsv file separated by comma",
)
@click.option(
    "--rna_tsvs",
    help="Path(s) to RNA-seq *_translating_ORFs.tsv file separated by comma",
)
@click.option(
    "--ribotricer_index",
    help=(
        "Path to the index file of ribotricer\n"
        "This file should be generated using ribotricer prepare-orfs (required for BAM input)"
    ),
)
@click.option("--prefix", help="Prefix to output file")
@click.option(
    "--filter_by_tx_annotation",
    help="transcript_type to filter regions by",
    type=str,
    default="protein_coding",
    show_default=True,
)
@click.option(
    "--phase_score_cutoff",
    type=float,
    default=CUTOFF,
    show_default=True,
    help="Phase score cutoff for determining active translation (required for BAM input)",
)
@click.option(
    "--min_valid_codons",
    type=int,
    default=MINIMUM_VALID_CODONS,
    show_default=True,
    help="Minimum number of codons with non-zero reads for determining active translation (required for BAM input)",
)
@click.option(
    "--sampling_ratio",
    type=float,
    default=0.33,
    show_default=True,
    help="Number of protein coding regions to sample per bootstrap",
)
@click.option(
    "--n_bootstraps",
    type=int,
    default=20000,
    show_default=True,
    help="Number of bootstraps",
)
def determine_cutoff_cmd(
    ribo_bams,
    rna_bams,
    ribo_tsvs,
    rna_tsvs,
    ribotricer_index,
    prefix,
    filter_by_tx_annotation,
    phase_score_cutoff,
    min_valid_codons,
    sampling_ratio,
    n_bootstraps,
):

    filter_by = _clean_input(filter_by_tx_annotation)

    ribo_stranded_protocols = []
    rna_stranded_protocols = []

    if ribo_bams and ribo_tsvs:
        sys.exit("Error: --ribo-bams and --rna_bams cannot be specified together")

    if rna_bams and rna_tsvs:
        sys.exit("Error: --rna-bams and --rna_tsvs cannot be specified together")

    if (ribo_bams and rna_tsvs) or (rna_bams and ribo_tsvs):
        sys.exit("Error: BAM and TSV inputs cannot be specified together")

    if ribotricer_index:
        if not os.path.isfile(ribotricer_index):
            sys.exit("Error: ribotricer index file not found")
    if ribo_bams:
        ribo_bams = _clean_input(ribo_bams)
        rna_bams = _clean_input(rna_bams)
    else:
        ribo_tsvs = _clean_input(ribo_tsvs)
        rna_tsvs = _clean_input(rna_tsvs)

    if ribo_bams and rna_bams:
        if not prefix:
            sys.exit("Error: --prefix required with BAM inputs")
        elif not ribotricer_index:
            sys.exit("Error: --ribotricer_index required with BAM inputs")
        else:
            determine_cutoff_bam(
                ribo_bams,
                rna_bams,
                ribotricer_index,
                prefix,
                ribo_stranded_protocols,
                rna_stranded_protocols,
                filter_by,
                sampling_ratio,
                n_bootstraps,
                phase_score_cutoff,
                min_valid_codons,
                report_all=True,
            )
    else:
        determine_cutoff_tsv(
            ribo_tsvs, rna_tsvs, filter_by, sampling_ratio, n_bootstraps
        )
