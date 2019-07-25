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


class Interval:
    """Class for interval
       All the intervals used in this project is 1-based and closed
    """

    def __init__(self, chrom=None, start=1, end=1, strand=None):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand

    def __eq__(self, other):
        """Override the default Equals behavior"""
        return (
            self.chrom == other.chrom
            and self.start == other.start
            and self.end == other.end
            and self.strand == other.strand
        )

    def __repr__(self):
        return "{}\t{}\t{}\t{}".format(self.chrom, self.start, self.end, self.strand)
