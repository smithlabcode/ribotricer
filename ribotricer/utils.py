"""Utilities for analysis"""

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

from collections import defaultdict
from typing import Final

import numpy as np
from numpy.typing import NDArray
from tqdm.autonotebook import tqdm

from .statistics import phasescore

tqdm.pandas()

CODON_TO_AA: Final[dict[str, str]] = {
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "ATG": "M",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AAC": "N",
    "AAT": "N",
    "AAA": "K",
    "AAG": "K",
    "AGC": "S",
    "AGT": "S",
    "AGA": "R",
    "AGG": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAC": "H",
    "CAT": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TTC": "F",
    "TTT": "F",
    "TTA": "L",
    "TTG": "L",
    "TAC": "Y",
    "TAT": "Y",
    "TAA": "-",
    "TAG": "-",
    "TGC": "C",
    "TGT": "C",
    "TGA": "-",
    "TGG": "W",
}


def parse_ccds(annotation: str, orfs: str, saveto: str) -> None:
    """Parse CCDS annotations.

    Parameters
    ----------
    annotation : str
        Path for annotation files of candidate ORFs.
    orfs : str
        Path for translating ORFs.
    saveto : str
        Output file name.
    """
    anno_oids = []
    real_oids = []
    ccds = defaultdict(list)
    with open(annotation) as anno:
        total_lines = len(["" for line in anno])
    with open(annotation) as anno:
        with tqdm(total=total_lines) as pbar:
            # Skip header
            anno.readline()
            for line in anno:
                pbar.update()
                if not line:
                    raise RuntimeError("annotation line cannot be empty")
                fields = line.split("\t")
                if len(fields) != 13:
                    raise RuntimeError(
                        "unexpected number of columns found for annotation file"
                    )
                oid = fields[0]
                gid = fields[4]
                ccds[gid].append(oid)
                anno_oids.append(oid)

    ccds_orfs = {}
    with open(orfs) as orf:
        total_lines = len(["" for line in orf])
    with open(orfs) as orf:
        with tqdm(total=total_lines) as pbar:
            # Skip header
            orf.readline()
            for line in orf:
                pbar.update()
                if not line:
                    raise RuntimeError("orf line cannot be empty")
                fields = line.split("\t")
                if len(fields) != 5:
                    raise RuntimeError(
                        "unexpected number of columns found for orf file"
                    )
                oid = fields[0]
                count = int(fields[2])
                corr = float(fields[3])
                pval = float(fields[4])
                ccds_orfs[oid] = (count, corr, pval)
                real_oids.append(oid)

    rename = dict(zip(anno_oids, real_oids))
    to_write = "Gene_ID\tCount\tPeriodicity\tPval\n"
    n_genes = 0
    for gid in ccds:
        n_genes += 1
        count, corr, pval = (0, 0, 1)
        for oid in ccds[gid]:
            oid = rename[oid]
            t_cnt, t_corr, t_pval = ccds_orfs[oid]
            if t_corr >= corr:
                count, corr, pval = (t_cnt, t_corr, t_pval)
        to_write += f"{gid}\t{count}\t{corr}\t{pval}\n"

    with open(saveto, "w") as output:
        output.write(to_write)


def benchmark(
    rna_file: str,
    ribo_file: str,
    prefix: str,
    cutoff: int = 5,
) -> None:
    """Benchmark RNA vs Ribo profiles.

    Parameters
    ----------
    rna_file : str
        Path to RNA profile file.
    ribo_file : str
        Path to Ribo profile file.
    prefix : str
        Output prefix.
    cutoff : int, optional
        Minimum coverage cutoff, by default 5.
    """
    rna: dict[str, list[int]] = {}
    ribo: dict[str, list[int]] = {}

    print("reading RNA profiles")
    with open(rna_file) as orf:
        total_lines = len(["" for line in orf])
    with open(rna_file) as orf:
        with tqdm(total=total_lines) as pbar:
            for line in orf:
                pbar.update()
                chrom, start, end, cat, gid, strand, cov = line.strip().split("\t")
                cov = [int(x) for x in cov.strip().split()]
                if strand == "-":
                    cov.reverse()
                ID = "_".join([chrom, start, end, cat, gid])
                rna[ID] = cov

    print("reading Ribo profiles")
    with open(ribo_file) as orf:
        total_lines = len(["" for line in orf])
    with open(ribo_file) as orf:
        with tqdm(total=total_lines) as pbar:
            for line in orf:
                pbar.update()
                chrom, start, end, cat, gid, strand, cov = line.strip().split("\t")
                cov = [int(x) for x in cov.strip().split()]
                if strand == "-":
                    cov.reverse()
                ID = "_".join([chrom, start, end, cat, gid])
                ribo[ID] = cov

    to_write = "ID\tribo_coh\trna_coh\tribo_cov\trna_cov\n"
    common_ids = set(ribo.keys()) & set(rna.keys())
    for ID in tqdm(common_ids):
        if sum(rna[ID]) >= cutoff and sum(ribo[ID]) >= cutoff and len(ribo[ID]) >= 10:
            rna_coh, rna_valid = phasescore(rna[ID])
            rna_cov = rna_valid / len(rna[ID])
            ribo_coh, ribo_valid = phasescore(ribo[ID])
            ribo_cov = ribo_valid / len(ribo[ID])

            to_write += f"{ID}\t{ribo_coh}\t{rna_coh}\t{ribo_cov}\t{rna_cov}\n"
    with open(f"{prefix}_results.txt", "w") as output:
        output.write(to_write)


def angle(cov: list[int], frame: int) -> tuple[list[float], int]:
    """Compute angles for coverage profile.

    Parameters
    ----------
    cov : list[int]
        Coverage profile.
    frame : int
        Frame offset.

    Returns
    -------
    tuple[list[float], int]
        Tuple of (angles, number of zero vectors).
    """
    ans: list[float] = []
    nzeros = 0
    cov = cov[frame:]
    i = 0
    while i + 2 < len(cov):
        if cov[i] == cov[i + 1] == cov[i + 2] == 0:
            i += 3
            nzeros += 1
        else:
            real = cov[i] - 0.5 * (cov[i + 1] + cov[i + 2])
            img = np.sqrt(3) / 2 * (cov[i + 1] - cov[i + 2])
            if real == img == 0:
                i += 3
            else:
                ans.append(np.arctan2(img, real))
                i += 3
    return ans, nzeros


def theta_dist(
    rna_file: str,
    ribo_file: str,
    frame_file: str,
    prefix: str,
    cutoff: int = 5,
) -> None:
    """Compute theta distribution from RNA and Ribo profiles.

    Parameters
    ----------
    rna_file : str
        Path to RNA profile file.
    ribo_file : str
        Path to Ribo profile file.
    frame_file : str
        Path to frame file.
    prefix : str
        Output prefix.
    cutoff : int, optional
        Minimum coverage cutoff, by default 5.
    """
    rna: dict[str, list[int]] = {}
    ribo: dict[str, list[int]] = {}
    frame: dict[str, int] = {}

    print("reading frame file")
    with open(frame_file) as frame_r:
        total_lines = len(["" for line in frame_r])
    with open(frame_file) as frame_r:
        with tqdm(total=total_lines) as pbar:
            for line in frame_r:
                pbar.update()
                name, frame_n, strand, length = line.strip().split("\t")
                frame[name] = int(frame_n)

    print("reading RNA profiles")
    with open(rna_file) as orf:
        total_lines = len(["" for line in orf])
    with open(rna_file) as orf:
        with tqdm(total=total_lines) as pbar:
            for line in orf:
                pbar.update()
                chrom, start, end, cat, gid, strand, cov = line.strip().split("\t")
                cov = [int(x) for x in cov.strip().split()]
                if strand == "-":
                    cov.reverse()
                ID = "_".join([chrom, start, end, cat, gid])
                rna[ID] = cov

    print("reading Ribo profiles")
    with open(ribo_file) as orf:
        total_lines = len(["" for line in orf])
    with open(ribo_file) as orf:
        with tqdm(total=total_lines) as pbar:
            for line in orf:
                pbar.update()
                chrom, start, end, cat, gid, strand, cov = line.strip().split("\t")
                cov = [int(x) for x in cov.strip().split()]
                if strand == "-":
                    cov.reverse()
                ID = "_".join([chrom, start, end, cat, gid])
                ribo[ID] = cov

    rna_angles = []
    ribo_angles = []
    rna_zeros = ribo_zeros = 0
    total_reads = 0
    total_length = 0
    total_ribo_reads = total_ribo_length = 0
    common_ids = set(ribo.keys()) & set(rna.keys())
    print("calculating angles")
    for ID in tqdm(common_ids):
        if sum(rna[ID]) >= cutoff and sum(ribo[ID]) >= cutoff and len(ribo[ID]) >= 10:
            cur_rna_angles, cur_rna_zeros = angle(rna[ID], frame[ID])
            rna_angles += cur_rna_angles
            rna_zeros += cur_rna_zeros
            cur_ribo_angles, cur_ribo_zeros = angle(ribo[ID], frame[ID])
            ribo_angles += cur_ribo_angles
            ribo_zeros += cur_ribo_zeros
            total_reads += sum(rna[ID])
            total_length += len(rna[ID])
            total_ribo_reads += sum(ribo[ID])
            total_ribo_length += len(ribo[ID])
    # generate theoretical theta dist from Poisson
    print("generating theoretical theta dist from Poisson")
    mean = total_reads / total_length
    poisson_cov = np.random.poisson(mean, total_length)
    poisson_angles, poisson_zeros = angle(poisson_cov, 0)
    with open(f"{prefix}_angle_stats.txt", "w") as output:
        output.write(f"total_rna_reads: {total_reads}\n")
        output.write(f"total_rna_ccds_length: {total_length}\n")
        output.write(f"total_ribo_reads: {total_ribo_reads}\n")
        output.write(f"total_ribo_ccds_length: {total_ribo_length}\n")
        output.write(f"mean reads: {mean}\n")
        output.write(f"rna zero vectors: {rna_zeros}\n")
        output.write(f"poisson zero vectors: {poisson_zeros}\n")
        output.write(f"ribo zero vectors: {ribo_zeros}\n")
    with open(f"{prefix}_rna_angles.txt", "w") as output:
        output.write("\n".join(map(str, rna_angles)))
    with open(f"{prefix}_ribo_angles.txt", "w") as output:
        output.write("\n".join(map(str, ribo_angles)))
    with open(f"{prefix}_poisson_angles.txt", "w") as output:
        output.write("\n".join(map(str, poisson_angles)))


def theta_rna(rna_file: str, prefix: str, cutoff: int = 10) -> None:
    """Compute theta distribution from RNA profiles.

    Parameters
    ----------
    rna_file : str
        Path to RNA profile file.
    prefix : str
        Output prefix.
    cutoff : int, optional
        Minimum coverage cutoff, by default 10.
    """
    rna: dict[str, list[int]] = {}

    print("reading RNA profiles")
    with open(rna_file) as orf:
        total_lines = len(["" for line in orf])
    with open(rna_file) as orf:
        with tqdm(total=total_lines) as pbar:
            # Skip header
            orf.readline()
            for line in orf:
                pbar.update()
                fields = line.split("\t")
                oid = fields[0]
                cov = fields[1]
                cov = cov[1:-1]
                cov = [int(x) for x in cov.split(", ")]
                if sum(cov) > cutoff:
                    rna[oid] = cov

    rna_angles = []
    for ID in tqdm(list(rna.keys())):
        rna_angles += angle(rna[ID], 0)
    with open(f"{prefix}_raw_rna_angles.txt", "w") as output:
        output.write("\n".join(map(str, rna_angles)))


def _nucleotide_to_codon_profile(
    profile: str | list[int],
) -> NDArray[np.int64]:
    """Summarize nucleotide profile to a codon level profile.

    Parameters
    ----------
    profile : str | list[int]
        Nucleotide profile as string or list.

    Returns
    -------
    NDArray[np.int64]
        Codon level profile.
    """
    if isinstance(profile, str):
        profile = eval(profile)
    profile_arr = np.array(profile)
    codon_profile: NDArray[np.int64] = np.add.reduceat(
        profile_arr, range(0, len(profile_arr), 3)
    )
    return codon_profile


def summarize_profile_to_codon_level(detected_orfs: str, saveto: str) -> None:
    """Collapse nucleotide level profiles in ribotricer to codon level.

    Parameters
    ----------
    detected_orfs : str
        Path to ribotricer detect-orfs output.
    saveto : str
        Path to write output to.
    """
    with open(saveto, "w") as fout:
        fout.write("ORF_ID\tcodon_profile\n")
        with open(detected_orfs) as fin:
            # Skip header
            fin.readline()
            for line in fin:
                fields = line.strip().split("\t")
                oid, otype, status = fields[:3]
                gene_id, gene_name, gene_type = fields[9:12]
                chrom, strand, start_codon, profile = fields[12:]

                profile_stripped = profile.strip()[1:-1].split(", ")
                if profile_stripped[0]:
                    profile = np.array(list(map(int, profile_stripped)))
                codon_profile = np.add.reduceat(profile, range(0, len(profile), 3))
                fout.write(f"{oid}\t{list(codon_profile)}\n")


def translate(seq: str) -> str:
    """Translate a given nucleotide sequence to an amino acid sequence.

    Parameters
    ----------
    seq : str
        Nucleotide sequence.

    Returns
    -------
    str
        Translated sequence of amino acids.
    """
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i : i + 3]
            protein += CODON_TO_AA[codon]
    return protein


def learn_ribotricer_cutoff(roc_input_file: str) -> tuple[float, float]:
    """Learn ribotricer phase score cutoff.

    Parameters
    ----------
    roc_input_file : str
        Path to ROC file generated using ribotricer benchmark.

    Returns
    -------
    tuple[float, float]
        Tuple of (cutoff, fscore).
    """
    import pandas as pd
    from sklearn.metrics import precision_recall_fscore_support

    data = pd.read_csv(roc_input_file, sep="\t")
    ribotricer_scores = data.ribotricer
    truth = data.truth
    precision_recall_fscore_support_list: list[list[float]] = []

    cutoffs = np.linspace(0, 1, 1000)

    for cutoff_val in cutoffs:
        predicted = np.where(ribotricer_scores > cutoff_val, 1, 0)
        s = precision_recall_fscore_support(
            truth, predicted, average="binary", pos_label=1
        )
        precision_recall_fscore_support_list.append([cutoff_val] + list(s))
    precision_recall_fscore_support_df = pd.DataFrame(
        precision_recall_fscore_support_list,
        columns=["cutoff", "precision", "recall", "fscore", "cutoff2"],
    )
    cutoff = precision_recall_fscore_support_df.loc[
        precision_recall_fscore_support_df["fscore"].idxmax()
    ].cutoff.iloc[0]
    fscore = precision_recall_fscore_support_df["fscore"].max()
    return cutoff, fscore
