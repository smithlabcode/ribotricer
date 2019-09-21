"""Command line interface for ribotricer
"""
# Part of ribotricer software
#
# Copyright (C) 2019 Wenzheng Li, Saket Choudhary and Andrew D Smith
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
from .const import CUTOFF
from .const import MINIMUM_VALID_CODONS

from .count_orfs import count_orfs
from .count_orfs import count_orfs_codon
from .detect_orfs import detect_orfs
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
        report_all,
    )


###################### count-orfs function #########################################
@cli.command(
    "count-orfs",
    context_settings=CONTEXT_SETTINGS,
    help="Count reads for detected ORFs",
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
    help="Count reads for detected ORFs",
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
@click.option("--saveto", help="Path to output file", required=True)
def orf_seq_cmd(ribotricer_index, fasta, saveto):
    if not os.path.isfile(ribotricer_index):
        sys.exit("Error: ribotricer index file not found")

    if not os.path.isfile(fasta):
        sys.exit("Error: fasta file not found")

    orf_seq(ribotricer_index, fasta, saveto)
