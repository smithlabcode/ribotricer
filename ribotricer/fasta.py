"""Process fasta files"""

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

import os
import warnings
from collections import OrderedDict
from typing import TYPE_CHECKING

from pyfaidx import Fasta

if TYPE_CHECKING:
    from .interval import Interval


class FastaReader:
    """Class for reading and querying fasta file.

    Attributes
    ----------
    fasta_location : str
        Path to the fasta file.
    fasta : Fasta
        pyfaidx Fasta object.
    """

    def __init__(self, fasta_location: str) -> None:
        """Initialize FastaReader.

        Parameters
        ----------
        fasta_location : str
            Path to fasta file.

        Raises
        ------
        Exception
            If the fasta file cannot be read.
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

    def query(self, intervals: list[Interval]) -> list[str]:
        """Query regions for sequence.

        Parameters
        ----------
        intervals : list[Interval]
            The intervals for fasta (one-based and full-closed).

        Returns
        -------
        list[str]
            An array containing sequences for each Interval.
            This function is agnostic of the strand information,
            the position in the scores corresponds to the interval.

        Raises
        ------
        Exception
            If start or end position exceeds chromosome length.
        """
        sequences: list[str] = []
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
                        "Chromosome start point exceeds chromosome length: {}>{}".format(
                            i.start, chrom_length
                        )
                    )
                elif i.end > chrom_length:
                    raise Exception(
                        "Chromosome end point exceeds chromosome length: {}>{}".format(
                            i.end, chrom_length
                        )
                    )
                seq = self.fasta.get_seq(i.chrom, i.start, i.end)
                sequences.append(seq)
        return sequences

    def complement(self, seq: str) -> str:
        """Complement a FASTA sequence.

        Parameters
        ----------
        seq : str
            String fasta sequence.

        Returns
        -------
        str
            Complement of input fasta.
        """
        complement_letters = {"A": "T", "C": "G", "T": "A", "G": "C"}
        seq = seq.upper()
        comp: list[str] = []
        for nuc in seq:
            if nuc in complement_letters:
                comp.append(complement_letters[nuc])
            else:
                comp.append(nuc)
        return "".join(comp)

    def reverse_complement(self, seq: str) -> str:
        """Reverse-complement a FASTA sequence.

        Parameters
        ----------
        seq : str
            String fasta sequence.

        Returns
        -------
        str
            Reverse complement of input fasta.
        """
        seq = seq.upper()
        return self.complement(seq)[::-1]

    @property
    def chromosomes(self) -> OrderedDict[str, int]:
        """Return list of chromosome and their sizes as in the fasta file.

        Returns
        -------
        OrderedDict[str, int]
            Dictionary with {"chr": length} format.
        """
        chroms: OrderedDict[str, int] = OrderedDict()
        for chrom in list(self.fasta.keys()):
            chroms[chrom] = len(self.fasta[chrom])
        return chroms
