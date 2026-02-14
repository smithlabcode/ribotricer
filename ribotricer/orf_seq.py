"""Generate sequences for ribotricer annotation"""

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
from typing import Final

import pandas as pd
from tqdm.autonotebook import tqdm

from .fasta import FastaReader
from .interval import Interval

tqdm.pandas()

# Codon to amino acid translation table
CODON_TABLE: Final[dict[str, str]] = {
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
    "TAA": "_",
    "TAG": "_",
    "TGC": "C",
    "TGT": "C",
    "TGA": "_",
    "TGG": "W",
}


def translate_nt_to_aa(seq: str) -> str:
    """Translate nucleotide sequence to amino acid sequence.

    Parameters
    ----------
    seq : str
        Nucleotide sequence.

    Returns
    -------
    str
        Amino acid sequence.
    """
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i : i + 3]
            if "N" in codon:
                protein += "X"
            elif codon not in CODON_TABLE:
                sys.stderr.write(
                    "Found unknown codon {}. Substituting with X..\n".format(codon)
                )
                protein += "X"
            else:
                protein += CODON_TABLE[codon]
    return protein


def orf_seq(
    ribotricer_index: str,
    genome_fasta: str,
    saveto: str,
    translate: bool = False,
) -> None:
    """Generate sequence for ribotricer annotation.

    Parameters
    ----------
    ribotricer_index : str
        Path to ribotricer generated annotation.
    genome_fasta : str
        Path to genome fasta.
    saveto : str
        Path to output file.
    translate : bool, optional
        Whether to translate to protein sequence, by default False.
    """
    fasta = FastaReader(genome_fasta)
    annotation_df = pd.read_csv(ribotricer_index, sep="\t")

    with open(saveto, "w") as fh:
        fh.write("ORF_ID\tsequence\n")
        for idx, row in tqdm(annotation_df.iterrows(), total=annotation_df.shape[0]):
            chrom = str(row.chrom)
            orf_id = row.ORF_ID
            coordinates = row.coordinate.split(",")
            strand = row.strand
            intervals: list[Interval] = []
            seq = ""

            for coordinate in coordinates:
                start, stop = coordinate.split("-")
                start_pos = int(start)
                stop_pos = int(stop)
                interval = Interval(chrom, start_pos, stop_pos, strand)
                intervals.append(interval)

            seq = ("").join(fasta.query(intervals))
            if strand == "-":
                seq = fasta.reverse_complement(seq)
            if translate:
                if len(seq) % 3 != 0:
                    sys.stderr.write(
                        f"WARNING: Sequence length with ORF ID '{orf_id}' is not "
                        "a multiple of three. Output sequence might be "
                        "truncated.\n"
                    )
                    seq = seq[0 : (len(seq) // 3) * 3]
                seq = translate_nt_to_aa(seq)
            fh.write("{}\t{}\n".format(orf_id, seq))
