"""Functions for finding all candidate ORFs"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import warnings

from collections import Counter
from collections import defaultdict

import pysam
from tqdm import *
import numpy as np
import pandas as pd

from .fasta import FastaReader
from .gtf import GTFReader
from .interval import Interval
from .common import merge_intervals
from .orf import ORF


def tracks_to_ivs(tracks):
    """
    Parameters
    ----------
    tracks: List[GTFTrack]
            list of gtf tracks
    
    Returns
    -------
    intervals: List[Interval]
               list of Interval
    """
    chrom = {track.chrom for track in tracks}
    strand = {track.strand for track in tracks}
    if len(chrom) != 1 or len(strand) != 1:
        print('fail to fetch seq: inconsistent chrom or strand')
        return None
    chrom = list(chrom)[0]
    strand = list(strand)[0]
    intervals = [
        Interval(chrom, track.start, track.end, strand) for track in tracks
    ]
    intervals = merge_intervals(intervals)
    return intervals


def transcript_to_genome_iv(start, end, intervals, reverse=False):
    """
    Parameters
    ----------
    start: int
           start position in transcript
    end: int
         end position in transcript
    intervals: List[Interval]
               coordinate in genome
    reverse: bool
             whether if it is on the reverse strand

    Returns
    -------
    ivs: List[Interval]
         the coordinate for start, end in genome
    """
    total_len = sum(i.end - i.start + 1 for i in intervals)
    if reverse:
        start, end = total_len - end - 1, total_len - start - 1
    ivs = []
    start_genome = None
    end_genome = None

    ### find start in genome
    cur = 0
    for i in intervals:
        i_len = i.end - i.start + 1
        if cur + i_len > start:
            start_genome = i.start + start - cur
            break
        cur += i_len

    ### find end in genome
    cur = 0
    for i in intervals:
        i_len = i.end - i.start + 1
        if cur + i_len > end:
            end_genome = i.start + end - cur
            break
        cur += i_len

    ### find overlap with (start_genome, end_genome)
    for i in intervals:
        s = max(i.start, start_genome)
        e = min(i.end, end_genome)
        if s <= e:
            ivs.append(Interval(i.chrom, s, e, i.strand))
    return ivs


def fetch_seq(fasta, tracks):
    """
    Parameters
    ----------
    fasta: FastaReader
           instance of FastaReader
    tracks: List[GTFTrack]
            list of gtf track

    Returns
    -------
    merged_seq: str
                combined seqeunce for the region
    """
    intervals = tracks_to_ivs(tracks)
    if not isinstance(fasta, FastaReader):
        fasta = FastaReader(fasta)
    sequences = fasta.query(intervals)
    merged_seq = ''.join(sequences)
    strand = tracks[0].strand
    if strand == '-':
        return fasta.reverse_complement(merged_seq)
    return merged_seq


def search_orfs(fasta, intervals):
    """
    Parameters
    ----------
    fasta: FastaReader
           instance of FastaReader
    intervals: List[Interval]
               list of intervals

    Returns
    -------
    orfs: list
          list of (List[Interval], seq, leader, trailer)
            list of intervals for candidate ORF
            seq: sequence for the candidate ORF
            leader: sequence upstream of the ORF
            trailer: sequence downstream of the ORF
    """
    if not intervals:
        return []

    orfs = []
    if not isinstance(fasta, FastaReader):
        fasta = FastaReader(fasta)
    intervals = merge_intervals(intervals)
    sequences = fasta.query(intervals)
    merged_seq = ''.join(sequences)
    reverse = False
    strand = intervals[0].strand
    if strand == '-':
        merged_seq = fasta.reverse_complement(merged_seq)
        reverse = True
    start_codons = set([
        'ATG', 'TTG', 'CTG', 'GTG', 'AAG', 'AGG', 'ACG', 'ACG', 'ATA', 'ATT',
        'ATC'
    ])
    stop_codons = set(['TAG', 'TAA', 'TGA'])
    for sc in start_codons:
        cur = 0
        while cur < len(merged_seq):
            start = merged_seq.find(sc, cur)
            if start == -1:
                break
            cur = start + 1
            for i in range(start, len(merged_seq), 3):
                if merged_seq[i:i + 3] in stop_codons:
                    ### found orf
                    ivs = transcript_to_genome_iv(start, i + 2, intervals,
                                                  reverse)
                    seq = merged_seq[start:i]
                    leader = merged_seq[:start]
                    trailer = merged_seq[i:]
                    if ivs:
                        orfs.append((ivs, seq, leader, trailer))
                    break
    return orfs


def check_orf_type(orf):
    return 'novel'


def prepare_orfs(gtf, fasta, prefix, min_orf_length, start_codons,
                 stop_codons):
    """
    Parameters
    ----------
    gtf: GTFReader
         instance of GTFReader
    fasta: FastaReader
           instance of FastaReader
    prefix: str
            prefix for output file
    min_orf_length: int
                    minimum length (nts) of ORF to include
    start_codons: set
                  set of start codons
    stop_codons: set
                 set of stop codons
    """

    candidate_orfs = []
    if not isinstance(gtf, GTFReader):
        gtf = GTFReader(gtf)
    if not isinstance(fasta, FastaReader):
        fasta = FastaReader(fasta)

    print('preparing candidate ORFs...')

    ### process CDS gtf
    print('searching cds...')
    cds_orfs = []
    for gid in tqdm(gtf.cds):
        for tid in gtf.cds[gid]:
            tracks = gtf.cds[gid][tid]
            orf = ORF.from_tracks(tracks, 'CDS')
            if orf:
                cds_orfs.append(orf)

    for tid in tqdm(gtf.transcript):
        tracks = gtf.transcript[tid]
        ttype = tracks[0].transcript_type
        gid = tracks[0].gene_id
        gname = tracks[0].gene_name
        gtype = tracks[0].gene_type
        chrom = tracks[0].chrom
        strand = tracks[0].strand
        ivs = tracks_to_ivs(tracks)
        orfs = search_orfs(fasta, ivs)
        for ivs, seq, leader, trailer in orfs:
            otype = check_orf_type(gid, tid, ivs)
            orf = ORF(otype, tid, ttype, gid, gname, gtype, chrom, strand, ivs,
                      seq, leader, trailer)
            candidate_orfs.append(orf)

    ### save to file
    print('saving candidate ORFs file...')
    to_write = ('ORF_ID\tORF_type\ttranscript_id\ttranscript_type'
                '\tgene_id\tgene_name\tgene_type\tchrom'
                '\tstrand\tcoordinate\tseq\tleader\ttrailer\n')
    formatter = '{}\t' * 12 + '{}\n'
    for orf in tqdm(cds_orfs):
        coordinate = ','.join(
            ['{}-{}'.format(iv.start, iv.end) for iv in orf.intervals])
        to_write += formatter.format(orf.oid, orf.category, orf.tid, orf.ttype,
                                     orf.gid, orf.gname, orf.gtype, orf.chrom,
                                     orf.strand, coordinate, orf.seq,
                                     orf.leader, orf.trailer)

    with open('{}_candidate_orfs.tsv'.format(prefix), 'w') as output:
        output.write(to_write)
