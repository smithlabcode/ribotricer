"""process fasta files"""
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

import os
import warnings

from pyfaidx import Fasta


class FastaReader:
    """Class for reading and querying fasta file."""

    def __init__(self, fasta_location):
        """
        Parameters
        ---------
        fasta_location : string
                         Path to fasta file

        """
        self.fasta_location = fasta_location
        try:
            self.fasta = Fasta(fasta_location, as_raw=True, sequence_always_upper=True)
        except Exception as e:
            raise Exception(
                "Error reading fasta file {} : {}".format(
                    os.path.abspath(self.fasta_location), e
                )
            )

    def query(self, intervals):
        """ Query regions for sequence.

        Parameters
        ----------
        intervals: list of Interval
                   The intervals for fasta is one-based and full-closed

        Returns
        -------
        sequences: list(str)
                   An array containing scores for each Interval
                   This function is agnostic of the strand information,
                   the position in the scores is corresponding to the interval

        .. currentmodule:: .FastaReader
        .. autosummary::
            .FastaReader

        """
        sequences = []
        chrom_lengths = self.chromosomes
        for i in intervals:
            if i.chrom not in list(chrom_lengths.keys()):
                warnings.warn(
                    "Chromosome {} does not appear in the fasta".format(i.chrom),
                    UserWarning,
                )
            else:
                chrom_length = chrom_lengths[i.chrom]
                if i.start > chrom_length:
                    raise Exception(
                        "Chromsome start point exceeds chromosome length: {}>{}".format(
                            i.start, chrom_length
                        )
                    )
                elif i.end > chrom_length:
                    raise Exception(
                        "Chromsome end point exceeds chromosome length: {}>{}".format(
                            i.end, chrom_length
                        )
                    )
                seq = self.fasta.get_seq(i.chrom, i.start, i.end)
                sequences.append(seq)
        return sequences

    def complement(self, seq):
        """Complement a FASTA sequence.

        Parameters
        ----------
        seq: str
            String fasta sequence


        Returns
        -------
        complement_seq: str
                        complemenet of input fasta
        """
        complement_letters = {"A": "T", "C": "G", "T": "A", "G": "C"}
        seq = seq.upper()
        comp = []
        for nuc in seq:
            if nuc in complement_letters:
                comp.append(complement_letters[nuc])
            else:
                comp.append(nuc)
        return "".join(comp)

    def reverse_complement(self, seq):
        """Reverse-complment a FASTA sequence.

        Parameters
        ----------
        seq: str
            String fasta sequence


        Returns
        -------
        complement_seq: str
                        complemenet of input fasta
        """
        seq = seq.upper()
        return self.complement(seq)[::-1]

    @property
    def chromosomes(self):
        """Return list of chromsome and their sizes
        as in the fasta file.

        Returns
        -------
        chroms : dict
                 Dictionary with {"chr": "Length"} format


        .. currentmodule:: .FastaReader
        .. autosummary::
            .FastaReader
        """
        chroms = {}
        for chrom in list(self.fasta.keys()):
            chroms[chrom] = len(self.fasta[chrom])
        return chroms
