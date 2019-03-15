"""Utilities for reading GTF file
"""

from collections import defaultdict
from tqdm import *


class GTFTrack:
    """Class for feature in GTF file."""

    standards = {
        'gene_biotype': 'gene_type',
        'transcript_biotype': 'transcript_type'
    }

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
                if k in GTFTrack.standards:
                    k = GTFTrack.standards[k]
                setattr(self, k, v.strip('"'))
        if not hasattr(self, 'gene_name') and hasattr(self, 'gene_id'):
            setattr(self, 'gene_name', self.gene_id)
        if not hasattr(self, 'transcript_name') and hasattr(
                self, 'transcript_id'):
            setattr(self, 'transcript_name', self.transcript_id)

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

        if feature not in ['exon', 'cds']:
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
        self.transcript = defaultdict(list)
        self.cds = defaultdict(lambda: defaultdict(list))
        # print('reading GTF file...')
        with open(self.gtf_location, 'r') as gtf:
            total_lines = len(['' for line in gtf])
        with open(self.gtf_location, 'r') as gtf:
            with tqdm(total=total_lines) as pbar:
                for line in gtf:
                    pbar.update()
                    track = GTFTrack.from_string(line)
                    if track is None:
                        continue
                    try:
                        gid = track.gene_id
                        tid = track.transcript_id
                    except AttributeError:
                        print('missing gene or transcript id {}:{}-{}'.format(
                            track.chrom, track.start, track.end))
                        continue

                    if track.feature == 'exon':
                        self.transcript[tid].append(track)
                    elif track.feature == 'cds':
                        self.cds[gid][tid].append(track)
