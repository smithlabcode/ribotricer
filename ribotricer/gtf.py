"""Utilities for reading GTF file"""

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
from typing import Any, ClassVar

from tqdm.autonotebook import tqdm

tqdm.pandas()


class GTFTrack:
    """Class for feature in GTF file.

    Attributes
    ----------
    chrom : str
        Chromosome name.
    source : str
        Source of the annotation.
    feature : str
        Feature type (e.g., 'exon', 'cds').
    start : int
        Start position (1-based).
    end : int
        End position (1-based).
    score : str
        Score field.
    strand : str
        Strand ('+' or '-').
    frame : str
        Reading frame.
    """

    standards: ClassVar[dict[str, str]] = {
        "gene_biotype": "gene_type",
        "transcript_biotype": "transcript_type",
    }

    def __init__(
        self,
        chrom: str,
        source: str,
        feature: str,
        start: int,
        end: int,
        score: str,
        strand: str,
        frame: str,
        attribute: str,
    ) -> None:
        self.chrom = chrom
        self.source = source
        self.feature = feature
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.frame = frame

        # Parse attributes
        for att in attribute.split(";"):
            if len(att.split()) == 2:
                k, v = att.strip().split()
                if k in GTFTrack.standards:
                    k = GTFTrack.standards[k]
                setattr(self, k, v.strip('"'))

        if not hasattr(self, "gene_name") and hasattr(self, "gene_id"):
            setattr(self, "gene_name", self.gene_id)
        if not hasattr(self, "transcript_name") and hasattr(self, "transcript_id"):
            setattr(self, "transcript_name", self.transcript_id)
        if not hasattr(self, "transcript_type") and not hasattr(
            self, GTFTrack.standards["transcript_biotype"]
        ):
            # transcript_type not set so set it to "assumed_protein_coding".
            setattr(self, "transcript_type", "assumed_protein_coding")
        if not hasattr(self, "gene_type") and hasattr(self, "transcript_type"):
            setattr(self, "gene_type", self.transcript_type)

    # Type hints for dynamically set attributes
    gene_id: str
    gene_name: str
    gene_type: str
    transcript_id: str
    transcript_name: str
    transcript_type: str

    @classmethod
    def from_string(cls, line: str) -> GTFTrack | None:
        """Parse a GTF line into a GTFTrack.

        Parameters
        ----------
        line : str
            One line in GTF file.

        Returns
        -------
        GTFTrack | None
            Parsed track object, or None if line should be skipped.

        Notes
        -----
        This method follows the fails-fast strategy and
        hence uses multiple returns, ultimately returning
        a line from the GTF parsed into a feature (chrom, start end etc.)
        """
        if line.startswith("#"):
            return None

        fields = line.strip().split("\t")
        if len(fields) != 9:
            print("mal-formatted GTF file")
            return None

        chrom = fields[0]
        source = fields[1]
        feature = fields[2].lower()
        start = int(fields[3])
        end = int(fields[4])
        score = fields[5]
        strand = fields[6]
        frame = fields[7]
        attribute = fields[8]

        if feature not in ["exon", "cds"]:
            return None

        return cls(chrom, source, feature, start, end, score, strand, frame, attribute)

    def __repr__(self) -> str:
        return str(self.__dict__)


class GTFReader:
    """Class for reading and parsing gtf file.

    Attributes
    ----------
    gtf_location : str
        Path to GTF file.
    transcript : defaultdict[str, list[GTFTrack]]
        Dictionary mapping transcript IDs to their exon tracks.
    cds : defaultdict[str, defaultdict[str, list[GTFTrack]]]
        Dictionary mapping gene IDs to transcript IDs to CDS tracks.
    """

    def __init__(self, gtf_location: str) -> None:
        """Initialize GTFReader.

        Parameters
        ----------
        gtf_location : str
            Path to GTF file.
        """
        self.gtf_location = gtf_location
        self.transcript: defaultdict[str, list[GTFTrack]] = defaultdict(list)
        self.cds: defaultdict[str, defaultdict[str, list[GTFTrack]]] = defaultdict(
            lambda: defaultdict(list)
        )

        with open(self.gtf_location, "r") as gtf:
            total_lines = len(["" for line in gtf])
        with open(self.gtf_location, "r") as gtf:
            with tqdm(total=total_lines, unit="lines", leave=False) as pbar:
                for line in gtf:
                    pbar.update()
                    track = GTFTrack.from_string(line)
                    if track is not None:
                        try:
                            gid = track.gene_id
                            tid = track.transcript_id
                        except AttributeError:
                            print(
                                "missing gene or transcript id {}:{}-{}".format(
                                    track.chrom, track.start, track.end
                                )
                            )
                        else:
                            if track.feature == "exon":
                                self.transcript[tid].append(track)
                            elif track.feature == "cds":
                                self.cds[gid][tid].append(track)
