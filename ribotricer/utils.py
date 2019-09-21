"""Utilities for analysis
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

from collections import defaultdict
import numpy as np
from tqdm import tqdm
from .statistics import coherence

CODON_TO_AA = {
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


def parse_ccds(annotation, orfs, saveto):
    """
    Parameters
    ----------
    annotation: str
                Path for annotation files of candidate ORFs
    orfs: str
          Path for translating ORFs
    saveto: str
          output file name
    """
    anno_oids = []
    real_oids = []
    ccds = defaultdict(list)
    with open(annotation, "r") as anno:
        total_lines = len(["" for line in anno])
    with open(annotation, "r") as anno:
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
    with open(orfs, "r") as orf:
        total_lines = len(["" for line in orf])
    with open(orfs, "r") as orf:
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

    rename = {x: y for (x, y) in zip(anno_oids, real_oids)}
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
        to_write += "{}\t{}\t{}\t{}\n".format(gid, count, corr, pval)

    with open(saveto, "w") as output:
        output.write(to_write)


def benchmark(rna_file, ribo_file, prefix, cutoff=5):

    rna = {}
    ribo = {}

    print("reading RNA profiles")
    with open(rna_file, "r") as orf:
        total_lines = len(["" for line in orf])
    with open(rna_file, "r") as orf:
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
    with open(ribo_file, "r") as orf:
        total_lines = len(["" for line in orf])
    with open(ribo_file, "r") as orf:
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
            rna_coh, rna_valid = coherence(rna[ID])
            rna_cov = rna_valid / len(rna[ID])
            ribo_coh, ribo_valid = coherence(ribo[ID])
            ribo_cov = ribo_valid / len(ribo[ID])

            to_write += "{}\t{}\t{}\t{}\t{}\n".format(
                ID, ribo_coh, rna_coh, ribo_cov, rna_cov
            )
    with open("{}_results.txt".format(prefix), "w") as output:
        output.write(to_write)


def angle(cov, frame):
    ans = []
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


def theta_dist(rna_file, ribo_file, frame_file, prefix, cutoff=5):

    rna = {}
    ribo = {}
    frame = {}

    print("reading frame file")
    with open(frame_file, "r") as frame_r:
        total_lines = len(["" for line in frame_r])
    with open(frame_file, "r") as frame_r:
        with tqdm(total=total_lines) as pbar:
            for line in frame_r:
                pbar.update()
                name, frame_n, strand, length = line.strip().split("\t")
                frame[name] = int(frame_n)

    print("reading RNA profiles")
    with open(rna_file, "r") as orf:
        total_lines = len(["" for line in orf])
    with open(rna_file, "r") as orf:
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
    with open(ribo_file, "r") as orf:
        total_lines = len(["" for line in orf])
    with open(ribo_file, "r") as orf:
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
    with open("{}_angle_stats.txt".format(prefix), "w") as output:
        output.write("total_rna_reads: {}\n".format(total_reads))
        output.write("total_rna_ccds_length: {}\n".format(total_length))
        output.write("total_ribo_reads: {}\n".format(total_ribo_reads))
        output.write("total_ribo_ccds_length: {}\n".format(total_ribo_length))
        output.write("mean reads: {}\n".format(mean))
        output.write("rna zero vectors: {}\n".format(rna_zeros))
        output.write("poisson zero vectors: {}\n".format(poisson_zeros))
        output.write("ribo zero vectors: {}\n".format(ribo_zeros))
    with open("{}_rna_angles.txt".format(prefix), "w") as output:
        output.write("\n".join(map(str, rna_angles)))
    with open("{}_ribo_angles.txt".format(prefix), "w") as output:
        output.write("\n".join(map(str, ribo_angles)))
    with open("{}_poisson_angles.txt".format(prefix), "w") as output:
        output.write("\n".join(map(str, poisson_angles)))


def theta_rna(rna_file, prefix, cutoff=10):

    rna = {}

    print("reading RNA profiles")
    with open(rna_file, "r") as orf:
        total_lines = len(["" for line in orf])
    with open(rna_file, "r") as orf:
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
    with open("{}_raw_rna_angles.txt".format(prefix), "w") as output:
        output.write("\n".join(map(str, rna_angles)))


def _nucleotide_to_codon_profile(profile):
    """Summarize nucleotid profile to a codon level profile"""
    if isinstance(profile, str):
        profile = eval(profile)
    profile = np.array(profile)
    codon_profile = np.add.reduceat(profile, range(0, len(profile), 3))
    return codon_profile


def summarize_profile_to_codon_level(detected_orfs, saveto):
    """Collapse nucleotide level profiles in ribotricer to codon leve.
  
  Parameters
  ----------
  ribotricer_output: string
                     Path to ribotricer detect-orfs output
  saveto: string
          Path to write output to
  """
    with open(saveto, "w") as fout:
        fout.write("ORF_ID\tcodon_profile\n")
        with open(detected_orfs, "r") as fin:
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
                fout.write("{}\t{}\n".format(oid, list(codon_profile)))


def translate(seq):
    """Translate a given nucleotide sequence to an amino acid sequence

  Parameters
  ----------
  seq: str
       Nucleotide seqeunce

  Returns
  -------
  protein: str
           Translated sequence of amino acids
  """

    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i : i + 3]
            protein += CODON_TO_AA[codon]
    return protein


def learn_ribotricer_cutoff(roc_input_file):
    """Learn ribotricer phase score cutoff 
  
  Parameters
  ----------
  roc_input_file: str
                  Path to ROC file generated using ribotricer benchmark

  Returns
  -------
  cutoff: float
          Recommended phase score cutoff  

  fscore: float
          Corresponding F1 score acheived at the determined cutoff
  """
    from sklearn.metrics import precision_recall_fscore_support
    import pandas as pd

    data = pd.read_csv(roc_input_file, sep="\t")
    ribotricer_scores = data.ribotricer
    truth = data.truth
    precision_recall_fscore_support_df = []

    cutoffs = np.linspace(0, 1, 1000)

    for cutoff in cutoffs:
        predicted = np.where(ribotricer_scores > cutoff, 1, 0)
        s = precision_recall_fscore_support(
            truth, predicted, average="binary", pos_label=1
        )
        precision_recall_fscore_support_df.append([cutoff] + list(s))
    precision_recall_fscore_support_df = pd.DataFrame(
        precision_recall_fscore_support_df,
        columns=["cutoff", "precision", "recall", "fscore", "cutoff"],
    )
    cutoff = precision_recall_fscore_support_df.loc[
        precision_recall_fscore_support_df["fscore"].idxmax()
    ].cutoff.iloc[0]
    fscore = precision_recall_fscore_support_df["fscore"].max()
    return cutoff, fscore
