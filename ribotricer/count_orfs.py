"""Utilities for translating ORF detection"""
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

from collections import defaultdict
from textwrap import wrap
from .orf import ORF

import numpy as np
import pandas as pd


def count_orfs(
    ribotricer_index, detected_orfs, features, outfile, report_all=False
):
    """
    Parameters
    ----------
    ribotricer_index: str
                      Path to the index file generated by ribotricer prepare_orfs
    detected_orfs: str
                   Path to the detected orfs file generated by ribotricer detect_orfs
    features: set
              set of ORF types, such as {annotated}
    prefix: str
            prefix for output file
    report_all: bool
                if True, all coverages will be exported
    """
    orf_index = {}
    read_counts = defaultdict(dict)
    with open(ribotricer_index, "r") as fin:
        # Skip header
        fin.readline()
        for line in fin:
            orf = ORF.from_string(line)
            if orf.category in features:
                orf_index[orf.oid] = orf
    with open(detected_orfs, "r") as fin:
        # Skip header
        fin.readline()
        for line in fin:
            fields = line.strip().split("\t")
            oid, otype, status = fields[:3]
            gene_id, gene_name, gene_type = fields[11:14]
            chrom, strand, start_codon, profile = fields[14:]
            if otype in features:
                # do not output 'nontranslating' events unless report_all is set
                if status != "nontranslating" or report_all:
                    intervals = orf_index[oid].intervals
                    coor = [
                        x
                        for iv in intervals
                        for x in range(iv.start, iv.end + 1)
                    ]
                    if strand == "-":
                        coor = coor[::-1]
                    profile_stripped = profile.strip()[1:-1].split(", ")
                    profile = list()
                    if profile_stripped[0]:
                        profile = list(map(int, profile_stripped))
                    for pos, cov in zip(coor, profile):
                        if pos not in read_counts[gene_id, gene_name]:
                            read_counts[gene_id, gene_name][pos] = cov

    # Output count table
    with open(outfile, "w") as fout:
        fout.write("gene_id\tcount\tlength\n")
        for gene_id, gene_name in sorted(read_counts):
            values = read_counts[gene_id, gene_name].values()
            length = len(values)
            total = sum(values)
            fout.write("{}\t{}\t{}\n".format(gene_id, total, length))


def count_orfs_codon(
    ribotricer_index,
    detected_orfs,
    features,
    ribotricer_index_fasta,
    prefix,
    report_all=False,
):
    """
    Create genewise codon summaries

    Parameters
    ----------
    ribotricer_index: str
                      Path to the index file generated by ribotricer prepare_orfs
    detected_orfs: str
                   Path to the detected orfs file generated by ribotricer detect_orfs
    features: set
              set of ORF types, such as {annotated}
    ribotricer_index_fasta: str
                            Path to fasta index generated using orf-seq
    prefix: str
           path to output file
    report_all: bool
                if True, all coverages will be exported
    """
    orf_index = {}
    fasta_df = pd.read_csv(ribotricer_index_fasta, sep="\t").set_index(
        "ORF_ID"
    )
    read_counts = defaultdict(dict)
    with open(ribotricer_index, "r") as fin:
        # Skip header
        fin.readline()
        for line in fin:
            orf = ORF.from_string(line)
            if orf.category in features:
                orf_index[orf.oid] = orf
    with open(detected_orfs, "r") as fin:
        # Skip header
        fin.readline()
        for line in fin:
            fields = line.strip().split("\t")
            oid, otype, status = fields[:3]
            gene_id, gene_name, gene_type = fields[11:14]
            chrom, strand, start_codon, profile = fields[14:]
            if otype in features:
                # do not output 'nontranslating' events unless report_all is set
                if status != "nontranslating" or report_all:
                    intervals = orf_index[oid].intervals
                    coor = [
                        x
                        for iv in intervals
                        for x in range(iv.start, iv.end + 1)
                    ]
                    codon_coor = [
                        x
                        for iv in intervals
                        for x in range(iv.start, iv.end + 1, 3)
                    ]
                    if strand == "-":
                        coor = coor[::-1]
                    profile_stripped = profile.strip()[1:-1].split(", ")
                    profile = list()
                    if profile_stripped[0]:
                        profile = list(map(int, profile_stripped))
                    # IMP: Skip profiles that are not 3n long to avoid errors
                    # downstream with sequenceu
                    if len(profile) % 3 != 0:
                        continue

                    codon_profile = np.add.reduceat(
                        profile, range(0, len(profile), 3)
                    ).tolist()
                    assert sum(codon_profile) == sum(profile)
                    codon_seq = str(fasta_df.loc[oid].sequence)
                    if not len(codon_seq) % 3 == 0:
                        print(oid, len(codon_seq))
                    codon_seq_partitioned = wrap(codon_seq, 3)
                    for pos, cov, codon_seq in zip(
                        codon_coor, codon_profile, codon_seq_partitioned
                    ):
                        if pos not in read_counts[gene_id, codon_seq]:
                            read_counts[gene_id, codon_seq][pos] = cov

    # Output count table
    with open("{}_genewise.tsv".format(prefix), "w") as fout:
        fout.write(
            "\t".join(
                "gene_id",
                "codon",
                "values",
                "mean_codon_coverage",
                "median_codon_coverage",
                "var_codon_coverage",
                "codon_occurences",
                "total_codon_coverage",
            )
            + "\n"
        )
        for gene_id, codon_seq in sorted(read_counts):
            values = list(read_counts[gene_id, codon_seq].values())
            codon_occurences = len(values)
            total_codon_coverage = sum(values)
            mean_codon_coverage = np.mean(values)
            median_codon_coverage = np.median(values)
            var_codon_coverage = np.var(values)
            fout.write(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    gene_id,
                    codon_seq,
                    values,
                    mean_codon_coverage,
                    median_codon_coverage,
                    var_codon_coverage,
                    codon_occurences,
                    total_codon_coverage,
                )
            )
    fout_df = pd.read_csv("{}_genewise.tsv".format(prefix), sep="\t")
    fout_df["per_codon_enrichment(total/n_occur)"] = (
        fout_df["total_codon_coverage"] / fout_df["codon_occurences"]
    )
    fout_df[
        "-log10_relative_enrichment(per_codon/total_gene_coverage)"
    ] = -np.log10(
        fout_df["per_codon_enrichment(total/n_occur)"]
        / fout_df.groupby("gene_id")["total_codon_coverage"].transform("sum")
    )
    # Overwrite
    fout_df.to_csv(
        "{}_genewise.tsv".format(prefix), sep="\t", index=False, header=True
    )
    # Remove infs
    fout_df = fout_df.replace([np.inf, -np.inf], np.nan)
    fout_df = fout_df.dropna()
    fout_df["relative_enrichment"] = fout_df[
        "per_codon_enrichment(total/n_occur)"
    ] / fout_df.groupby("gene_id")["total_codon_coverage"].transform("sum")

    relative_enrichment_median = pd.DataFrame(
        fout_df.groupby("codon")["relative_enrichment"].median()
    )
    relative_enrichment_median.columns = ["median_relative_enrichment"]

    relative_enrichment_mean = pd.DataFrame(
        fout_df.groupby("codon")["relative_enrichment"].mean()
    )
    relative_enrichment_mean.columns = ["mean_relative_enrichment"]

    relative_enrichment_var = pd.DataFrame(
        fout_df.groupby("codon")["relative_enrichment"].var()
    )
    relative_enrichment_var.columns = ["var_relative_enrichment"]

    relative_enrichment = relative_enrichment_mean.join(
        relative_enrichment_median
    ).join(relative_enrichment_var)
    relative_enrichment.index.name = "codon"

    relative_enrichment = relative_enrichment.reset_index()
    relative_enrichment.to_csv(
        "{}_codonwise.tsv".format(prefix), sep="\t", index=False, header=True
    )
