"""Utilities for translating learning phase-score cutoffs"""

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

from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from .common import mkdir_p, parent_dir
from .const import (
    CUTOFF,
    MINIMUM_DENSITY_OVER_ORF,
    MINIMUM_READS_PER_CODON,
    MINIMUM_VALID_CODONS,
    MINIMUM_VALID_CODONS_RATIO,
)
from .detect_orfs import detect_orfs


def determine_cutoff_tsv(
    ribo_tsvs: list[str],
    rna_tsvs: list[str],
    filter_by: list[str] | None = None,
    sampling_ratio: float = 0.33,
    reps: int = 10000,
) -> None:
    """Learn cutoff empirically from ribotricer generated ORF TSVs.

    Parameters
    ----------
    ribo_tsvs : list[str]
        List of filepaths of ribotricer generated *translating_ORFs.tsv
        for Ribo-seq samples.
    rna_tsvs : list[str]
        List of filepaths of ribotricer generated *translating_ORFs.tsv
        for RNA-seq samples.
    filter_by : list[str] | None, optional
        Transcript types to filter by, by default ['protein_coding'].
    sampling_ratio : float, optional
        Sampling ratio, by default 0.33.
    reps : int, optional
        Number of replicates, by default 10000.
    """
    if filter_by is None:
        filter_by = ["protein_coding"]
    ribo_df = pd.DataFrame()
    for tsv in ribo_tsvs:
        df = pd.read_csv(
            tsv,
            sep="\t",
            usecols=["ORF_ID", "ORF_type", "phase_score", "transcript_type"],
        )
        ribo_df = pd.concat([ribo_df, df])

    rna_df = pd.DataFrame()
    for tsv in rna_tsvs:
        df = pd.read_csv(
            tsv,
            sep="\t",
            usecols=["ORF_ID", "ORF_type", "phase_score", "transcript_type"],
        )
        rna_df = pd.concat([rna_df, df])

    filter_by = list(map(lambda x: x.lower(), filter_by))

    ribo_df_filtered = ribo_df.loc[ribo_df.ORF_type == "annotated"]
    ribo_df_filtered = ribo_df_filtered.loc[
        ribo_df_filtered.transcript_type.str.lower().isin(filter_by)
    ]

    rna_df_filtered = rna_df.loc[rna_df.ORF_type == "annotated"]
    rna_df_filtered = rna_df_filtered.loc[
        rna_df_filtered.transcript_type.str.lower().isin(filter_by)
    ]

    n_total_ribo = ribo_df_filtered.shape[0]
    n_total_rna = rna_df_filtered.shape[0]

    n_select_ribo = int(sampling_ratio * n_total_ribo)
    n_select_rna = int(sampling_ratio * n_total_rna)

    np.random.seed(42)
    ribo_indices = np.random.choice(range(n_total_ribo), (n_select_ribo, reps))
    rna_indices = np.random.choice(range(n_total_rna), (n_select_rna, reps))

    ribo_phase_scores_all = np.array(ribo_df_filtered.phase_score.tolist())
    rna_phase_scores_all = np.array(rna_df_filtered.phase_score.tolist())

    ribo_phase_median_samples = np.median(ribo_phase_scores_all[ribo_indices], axis=0)
    rna_phase_median_samples = np.median(rna_phase_scores_all[rna_indices], axis=0)
    diff_phase_median_samples = ribo_phase_median_samples - rna_phase_median_samples

    ribo_phase_score_mean = np.mean(ribo_phase_median_samples)
    ribo_phase_score_median = np.median(ribo_phase_median_samples)
    ribo_phase_score_sd = np.std(ribo_phase_median_samples)

    rna_phase_score_mean = np.mean(rna_phase_median_samples)
    rna_phase_score_median = np.median(rna_phase_median_samples)
    rna_phase_score_sd = np.std(rna_phase_median_samples)

    diff_phase_score_mean = np.mean(diff_phase_median_samples)
    diff_phase_score_median = np.median(diff_phase_median_samples)
    diff_phase_score_sd = np.std(diff_phase_median_samples)

    diff_all = ribo_phase_scores_all - rna_phase_scores_all
    diff_all_median = np.median(diff_all)
    diff_all_mean = np.mean(diff_all)
    diff_all_std = np.std(diff_all)

    print("sampling_ratio: {}".format(sampling_ratio))
    print("n_samples: {}".format(reps))

    print("ribo_phase_score_mean: {:.3f}".format(ribo_phase_score_mean))
    print("ribo_phase_score_median: {:.3f}".format(ribo_phase_score_median))
    print("ribo_phase_score_sd: {:.3f}".format(ribo_phase_score_sd))

    print("rna_phase_score_mean: {:.3f}".format(rna_phase_score_mean))
    print("rna_phase_score_median: {:.3f}".format(rna_phase_score_median))
    print("rna_phase_score_sd: {:.3f}".format(rna_phase_score_sd))

    print("diff_phase_score_sampled_mean: {:.3f}".format(diff_phase_score_mean))
    print("diff_phase_score_sampled_median: {:.3f}".format(diff_phase_score_median))
    print("diff_phase_score_sampled_sd: {:.3f}".format(diff_phase_score_sd))

    print("diff_phase_score_all_mean: {:.3f}".format(diff_all_mean))
    print("diff_phase_score_all_median: {:.3f}".format(diff_all_median))
    print("diff_phase_score_all_sd: {:.3f}".format(diff_all_std))

    print("recommended_cutoff: {:.3f}".format(diff_phase_score_median))


def determine_cutoff_bam(
    ribo_bams: list[str],
    rna_bams: list[str],
    ribotricer_index: str,
    prefix: str,
    ribo_stranded_protocols: list[str | None] | None = None,
    rna_stranded_protocols: list[str | None] | None = None,
    filter_by: list[str] | None = None,
    sampling_ratio: float = 0.33,
    reps: int = 10000,
    phase_score_cutoff: float = CUTOFF,
    min_valid_codons: int = MINIMUM_VALID_CODONS,
    report_all: bool = True,
) -> None:
    """Learn cutoff empirically from the given data.

    This uses the following steps:

    1. Run ribotricer using a cutoff of 0 for both RNA and Ribo samples
    2. For each output of ribotricer, find the median difference between RNA
       and Ribo-seq phase scores using the protein_coding annotated regions
       in the output.

    Parameters
    ----------
    ribo_bams : list[str]
        List of filepaths to Ribo-seq BAMs.
    rna_bams : list[str]
        List of filepaths to RNA-seq BAMs.
    ribotricer_index : str
        Path to ribotricer index file.
    prefix : str
        Output prefix.
    ribo_stranded_protocols : list[str | None] | None, optional
        List of 'yes/no/reverse', by default None.
    rna_stranded_protocols : list[str | None] | None, optional
        List of 'yes/no/reverse', by default None.
    filter_by : list[str] | None, optional
        Transcript types to filter by, by default ['protein_coding'].
    sampling_ratio : float, optional
        Sampling ratio, by default 0.33.
    reps : int, optional
        Number of replicates, by default 10000.
    phase_score_cutoff : float, optional
        Phase score cutoff, by default CUTOFF.
    min_valid_codons : int, optional
        Minimum valid codons, by default MINIMUM_VALID_CODONS.
    report_all : bool, optional
        Whether to report all ORFs, by default True.
    """
    if ribo_stranded_protocols is None:
        ribo_stranded_protocols = []
    if rna_stranded_protocols is None:
        rna_stranded_protocols = []
    if filter_by is None:
        filter_by = ["protein_coding"]
    if len(ribo_stranded_protocols) > 1:
        if len(ribo_stranded_protocols) != len(ribo_bams):
            sys.exit("Error: Ribo-seq protocol and bam file length mismatch")
    else:
        ribo_stranded_protocols = [None] * len(ribo_bams)
    if len(rna_stranded_protocols) > 1:
        if len(ribo_stranded_protocols) != len(ribo_bams):
            sys.exit("Error: Ribo-seq protocol and bam file length mismatch")
    else:
        rna_stranded_protocols = [None] * len(rna_bams)
        ribo_stranded_protocols = [None] * len(ribo_bams)

    sample_str = "samples"
    if len(rna_bams) == 1:
        sample_str = "sample"
    print(
        "Running ribotricer on {} Ribo-seq {} ..... \n".format(
            len(rna_bams), sample_str
        )
    )
    ribo_bams_renamed = dict(
        zip(ribo_bams, ["ribo_bam_{}".format(i + 1) for i in range(len(ribo_bams))])
    )
    rna_bams_renamed = dict(
        zip(rna_bams, ["rna_bam_{}".format(i + 1) for i in range(len(rna_bams))])
    )

    rna_tsvs = []
    ribo_tsvs = []
    for bam, stranded in zip(ribo_bams, ribo_stranded_protocols):
        bam_prefix = "{}__{}".format(prefix, ribo_bams_renamed[bam])
        mkdir_p(parent_dir(bam_prefix))
        detect_orfs(
            bam,
            ribotricer_index,
            bam_prefix,
            stranded,
            read_lengths=None,
            psite_offsets=None,
            phase_score_cutoff=0.0,
            min_valid_codons=MINIMUM_VALID_CODONS,
            min_reads_per_codon=MINIMUM_READS_PER_CODON,
            min_valid_codons_ratio=MINIMUM_VALID_CODONS_RATIO,
            min_density_over_orf=MINIMUM_DENSITY_OVER_ORF,
            report_all=report_all,
        )
        ribo_tsvs.append("{}_translating_ORFs.tsv".format(bam_prefix))
    print(
        "Running ribotricer on {} RNA-seq {} ..... \n".format(len(rna_bams), sample_str)
    )
    for bam, stranded in zip(rna_bams, rna_stranded_protocols):
        bam_prefix = "{}__{}".format(prefix, rna_bams_renamed[bam])
        mkdir_p(parent_dir(bam_prefix))
        detect_orfs(
            bam,
            ribotricer_index,
            bam_prefix,
            stranded,
            read_lengths=None,
            psite_offsets=None,
            phase_score_cutoff=0.0,
            min_valid_codons=MINIMUM_VALID_CODONS,
            min_reads_per_codon=MINIMUM_READS_PER_CODON,
            min_valid_codons_ratio=MINIMUM_VALID_CODONS_RATIO,
            min_density_over_orf=MINIMUM_DENSITY_OVER_ORF,
            report_all=report_all,
        )
        rna_tsvs.append("{}_translating_ORFs.tsv".format(bam_prefix))
    determine_cutoff_tsv(ribo_tsvs, rna_tsvs, filter_by, sampling_ratio, reps)
