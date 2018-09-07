"""Utilities for reading GTF file
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from collections import defaultdict
from bx.intervals.intersection import IntervalTree


class GTFTrack:
    """Class for feature in GTF file."""

    def __init__(self, seqname, source, feature, start, end, score, strand,
                 frame, attribute):
        self.seqname = seqname
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
                setattr(self, k, v)

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

        fields = line.strip().split()
        if len(fields) < 9:
            print('mal-formatted GTF file')
            return None

        seqname = fields[0]
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

        return cls(seqname, source, feature, start, end, score, strand, frame,
                   attribute)


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
        self.utr = defaultdict(list)
        with open(self.gtf_location, 'r') as gtf:
            for line in gtf:
                track = GTFTrack.from_string(line)
                if track is None:
                    continue
                if track.feature == 'gene':
                    self.gene[track.seqname].insert(track.start, track.end,
                                                    track.strand)
                    continue
                try:
                    gid = track.gene_id
                    tid = track.transcript_id
                except AttributeError:
                    print('missing gene or transcript id {}:{}-{}'.format(
                        track.seqname, track.start, track.end))
                    continue

                if track.feature == 'cds':
                    self.cds[track.gene_id][track.transcript_id].append(track)
                elif track.feature == 'utr':
                    self.utr[track.transcript_id].append(track)
