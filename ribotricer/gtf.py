"""Utilities for reading GTF file
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
from tqdm import tqdm


class GTFTrack(object):
    """Class for feature in GTF file."""

    standards = {"gene_biotype": "gene_type", "transcript_biotype": "transcript_type"}

    def __init__(
        self, chrom, source, feature, start, end, score, strand, frame, attribute
    ):
        self.chrom = chrom
        self.source = source
        self.feature = feature
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.frame = frame
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

    @classmethod
    def from_string(cls, line):
        """
        Parameters
        ----------
        line: string
              one line in gtf file

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

    def __repr__(self):
        return str(self.__dict__)


class GTFReader(object):
    """Class for reading and parseing gtf file."""

    def __init__(self, gtf_location):
        """
        Parameters
        ---------
        gtf_location : string
                       Path to gtf file
        """
        self.gtf_location = gtf_location
        self.transcript = defaultdict(list)
        self.cds = defaultdict(lambda: defaultdict(list))
        # print('reading GTF file...')
        with open(self.gtf_location, "r") as gtf:
            total_lines = len(["" for line in gtf])
        with open(self.gtf_location, "r") as gtf:
            with tqdm(total=total_lines) as pbar:
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
