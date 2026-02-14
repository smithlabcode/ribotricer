"""Utility for handling chromosome intervals"""

# Part of ribotricer software
#
# Copyright (C) 2020-2026 Saket Choudhary, Wenzheng Li, and Andrew D Smith
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


class Interval:
    """Class for interval.

    All the intervals used in this project are 1-based and closed.

    Attributes
    ----------
    chrom : str | None
        Chromosome name.
    start : int
        Start position (1-based, inclusive).
    end : int
        End position (1-based, inclusive).
    strand : str
        Strand ('+' or '-').
    """

    __slots__ = ("chrom", "start", "end", "strand")

    def __init__(
        self,
        chrom: str | None = None,
        start: int = 1,
        end: int = 1,
        strand: str = "+",
    ) -> None:
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand

    def __eq__(self, other: object) -> bool:
        """Override the default Equals behavior."""
        if not isinstance(other, Interval):
            return NotImplemented
        return (
            self.chrom == other.chrom
            and self.start == other.start
            and self.end == other.end
            and self.strand == other.strand
        )

    def __repr__(self) -> str:
        return f"{self.chrom}\t{self.start}\t{self.end}\t{self.strand}"

    def __len__(self) -> int:
        """Return the length of the interval."""
        return self.end - self.start + 1

    def __hash__(self) -> int:
        """Make Interval hashable."""
        return hash((self.chrom, self.start, self.end, self.strand))
