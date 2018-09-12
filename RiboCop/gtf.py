"""Utilities for reading GTF file
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from collections import defaultdict
from bx.intervals.intersection import IntervalTree
from tqdm import *


class GTFTrack:
    """Class for feature in GTF file."""

    def __init__(self, chrom, source, feature, start, end, score, strand,
                 frame, attribute):
        self.chrom = chrom
        self.source = source
        self.feature = feature
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.frame = frame
        for att in attribute.split(';'):
            if len(att.split()) == 2:
                k, v = att.strip().split()
                setattr(self, k, v.strip('"'))

    @classmethod
    def from_string(cls, line):
        """
        Parameters
        ----------
        line: string
              one line in gtf file
        """
        if line.startswith('#'):
            return None

        fields = line.strip().split('\t')
        if len(fields) != 9:
            print('mal-formatted GTF file')
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

        if feature not in ['gene', 'cds', 'utr']:
            return None

        return cls(chrom, source, feature, start, end, score, strand, frame,
                   attribute)

    def __repr__(self):
        return str(self.__dict__)


class GTFReader:
    """Class for reading and parseing gtf file."""

    def __init__(self, gtf_location):
        """
        Parameters
        ---------
        gtf_location : string
                       Path to gtf file

        """
        self.gtf_location = gtf_location
        self.gene = defaultdict(IntervalTree)
        self.cds = defaultdict(lambda: defaultdict(list))
        self.utr = defaultdict(lambda: defaultdict(list))
        with open(self.gtf_location, 'r') as gtf:
            for line in tqdm(gtf):
                track = GTFTrack.from_string(line)
                if track is None:
                    continue
                if track.feature == 'gene':
                    self.gene[track.chrom].insert(track.start, track.end,
                                                  track.strand)
                    continue
                try:
                    gid = track.gene_id
                    tid = track.transcript_id
                except AttributeError:
                    print('missing gene or transcript id {}:{}-{}'.format(
                        track.chrom, track.start, track.end))
                    continue

                if track.feature == 'cds':
                    self.cds[gid][tid].append(track)
                elif track.feature == 'utr':
                    self.utr[gid][tid].append(track)
