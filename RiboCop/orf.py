"""Utilities for translating ORF detection
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import warnings

from collections import Counter
from collections import defaultdict

import pysam
from tqdm import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

from .fasta import FastaReader
from .gtf import GTFReader
from .interval import Interval
from .common import is_read_uniq_mapping
from .common import merge_intervals
from .infer_protocol import infer_protocol


class PutativeORF:
    """Class for putative ORF."""

    def __init__(self,
                 category,
                 transcript_id,
                 transcript_type,
                 gene_id,
                 gene_name,
                 gene_type,
                 chrom,
                 strand,
                 intervals,
                 seq='',
                 leader='',
                 trailer=''):
        self.category = category
        self.tid = transcript_id
        self.ttype = transcript_type
        self.gid = gene_id
        self.gname = gene_name
        self.gtype = gene_type
        self.chrom = chrom
        self.strand = strand
        self.intervals = sorted(intervals, key=lambda x: x.start)
        start = intervals[0].start
        end = intervals[-1].end
        self.seq = seq
        self.oid = '{}_{}_{}_{}'.format(transcript_id, start, end, len(seq))
        self.leader = leader
        self.trailer = trailer

    @property
    def start_codon(self):
        if len(self.seq) < 3:
            return None
        return self.seq[:3]

    @classmethod
    def from_string(cls, line):
        """
        Parameters
        ----------
        line: string
              line for annotation file generated by prepare_orfs
        """
        if not line:
            print('annotation line cannot be empty')
            return None
        fields = line.strip().split('\t')
        if len(fields) != 13:
            print('unexpected number of columns found for annotation file')
            return None
        oid = fields[0]
        category = fields[1]
        tid = fields[2]
        ttype = fields[3]
        gid = fields[4]
        gname = fields[5]
        gtype = fields[6]
        chrom = fields[7]
        strand = fields[8]
        coordinate = fields[9]
        intervals = []
        for group in coordinate.strip().split(','):
            start, end = group.strip().split('-')
            start = int(start)
            end = int(end)
            intervals.append(Interval(chrom, start, end, strand))
        seq = fields[10]
        leader = fields[11]
        trailer = fields[12]
        return cls(category, tid, ttype, gid, gname, gtype, chrom, strand,
                   intervals)

    @classmethod
    def from_tracks(cls, tracks, category, seq='', leader='', trailer=''):
        """
        Parameters
        ----------
        tracks: list of GTFTrack
        """
        if not tracks:
            return None
        intervals = []
        tid = set()
        ttype = set()
        gid = set()
        gname = set()
        gtype = set()
        chrom = set()
        strand = set()
        for track in tracks:
            try:
                tid.add(track.transcript_id)
                ttype.add(track.transcript_type)
                gid.add(track.gene_id)
                gname.add(track.gene_name)
                gtype.add(track.gene_type)
                chrom.add(track.chrom)
                strand.add(track.strand)
                intervals.append(
                    Interval(track.chrom, track.start, track.end,
                             track.strand))
            except AttributeError:
                print('missing attribute {}:{}-{}'.format(
                    track.chrom, track.start, track.end))
                return None
        if (len(tid) != 1 or len(ttype) != 1 or len(gid) != 1
                or len(gname) != 1 or len(gtype) != 1 or len(chrom) != 1
                or len(strand) != 1):
            print('inconsistent tracks for one ORF')
            return None
        tid = list(tid)[0]
        ttype = list(ttype)[0]
        gid = list(gid)[0]
        gname = list(gname)[0]
        gtype = list(gtype)[0]
        chrom = list(chrom)[0]
        strand = list(strand)[0]
        return cls(category, tid, ttype, gid, gname, gtype, chrom, strand,
                   intervals, seq, leader, trailer)


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
            list of intervals for putative ORF
            seq: sequence for the putative ORF
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


def prepare_orfs(gtf, fasta, prefix):
    """
    Parameters
    ----------
    gtf: GTFReader
         instance of GTFReader
    fasta: FastaReader
           instance of FastaReader
    prefix: str
            prefix for output file

    Returns
    -------
    cds: List[PutativeORF]
         list of CDS
    uorfs: List[PutativeORF]
           list of upstream ORFs
    dorfs: List[PutativeORF]
           list of downstream ORFs
    """

    if not isinstance(gtf, GTFReader):
        gtf = GTFReader(gtf)
    if not isinstance(fasta, FastaReader):
        fasta = FastaReader(fasta)

    print('preparing putative ORFs...')

    ### process CDS gtf
    print('searching cds...')
    cds_orfs = []
    for gid in tqdm(gtf.cds):
        for tid in gtf.cds[gid]:
            tracks = gtf.cds[gid][tid]
            # seq = fetch_seq(fasta, tracks)
            orf = PutativeORF.from_tracks(tracks, 'CDS')
            if orf:
                cds_orfs.append(orf)

    ### process UTR gtf
    utr5 = defaultdict(list)
    utr3 = defaultdict(list)
    for gid in gtf.utr:
        ### find first cds and last cds for gene
        gene_cds = []
        for tid in gtf.cds[gid]:
            gene_cds += gtf.cds[gid][tid]
        if not gene_cds:
            print('fail to find CDS for UTR')
            continue
        first_cds = gene_cds[0]
        for gc in gene_cds:
            if gc.start < first_cds.start:
                first_cds = gc
        last_cds = gene_cds[-1]
        for gc in gene_cds:
            if gc.end > last_cds.end:
                last_cds = gc

        for tid in gtf.utr[gid]:
            for track in gtf.utr[gid][tid]:
                if track.start < first_cds.start:
                    if track.end >= first_cds.start:
                        track.end = first_cds.start - 1
                    if track.strand == '+':
                        utr5[tid].append(track)
                    else:
                        utr3[tid].append(track)
                elif track.end > last_cds.end:
                    if track.start <= last_cds.end:
                        track.start = last_cds.end + 1
                    if track.strand == '+':
                        utr3[tid].append(track)
                    else:
                        utr5[tid].append(track)

    uorfs = []
    print('searching uORFs...')
    for tid in tqdm(utr5):
        tracks = utr5[tid]
        ttype = tracks[0].transcript_type
        gid = tracks[0].gene_id
        gname = tracks[0].gene_name
        gtype = tracks[0].gene_type
        chrom = tracks[0].chrom
        strand = tracks[0].strand

        ivs = tracks_to_ivs(tracks)
        orfs = search_orfs(fasta, ivs)
        for ivs, seq, leader, trailer in orfs:
            orf = PutativeORF('uORF', tid, ttype, gid, gname, gtype, chrom,
                              strand, ivs, seq, leader, trailer)
            uorfs.append(orf)

    dorfs = []
    print('searching dORFs...')
    for tid in tqdm(utr3):
        tracks = utr3[tid]
        ttype = tracks[0].transcript_type
        gid = tracks[0].gene_id
        gname = tracks[0].gene_name
        gtype = tracks[0].gene_type
        chrom = tracks[0].chrom
        strand = tracks[0].strand

        ivs = tracks_to_ivs(tracks)
        orfs = search_orfs(fasta, ivs)
        for ivs, seq, leader, trailer in orfs:
            orf = PutativeORF('dORF', tid, ttype, gid, gname, gtype, chrom,
                              strand, ivs, seq, leader, trailer)
            dorfs.append(orf)

    ### save to file
    print('saving putative ORFs file...')
    to_write = ('ORF_ID\tORF_type\ttranscript_id\ttranscript_type'
                '\tgene_id\tgene_name\tgene_type\tchrom'
                '\tstrand\tcoordinate\tseq\tleader\ttrailer\n')
    formatter = '{}\t' * 12 + '{}\n'
    for orf in tqdm(cds_orfs + uorfs + dorfs):
        coordinate = ','.join(
            ['{}-{}'.format(iv.start, iv.end) for iv in orf.intervals])
        to_write += formatter.format(orf.oid, orf.category, orf.tid, orf.ttype,
                                     orf.gid, orf.gname, orf.gtype, orf.chrom,
                                     orf.strand, coordinate, orf.seq,
                                     orf.leader, orf.trailer)

    with open('{}_putative_orfs.tsv'.format(prefix), 'w') as output:
        output.write(to_write)

    return (cds_orfs, uorfs, dorfs)


def split_bam(bam, protocol, prefix, countby='5prime'):
    """Split bam by read length and strand

    Parameters
    ----------
    bam : str
          Path to bam file
    protocol: str
          Experiment protocol [forward, reverse]
    prefix: str
            prefix for output files
    countby: str
             5prime or 3prime to count the read

    Returns
    -------
    alignments: dict(dict(Counter))
                bam split by length, strand, (chrom, pos)
    read_lengths: dict
                  key is the length, value is the number of reads
    """
    alignments = defaultdict(lambda: defaultdict(Counter))
    read_lengths = defaultdict(int)
    qcfail = duplicate = secondary = unmapped = multi = valid = 0
    bam = pysam.AlignmentFile(bam, 'rb')
    total_count = bam.count()
    with tqdm(total=total_count) as pbar:
        for r in bam.fetch(until_eof=True):

            if r.is_qcfail:
                qcfail += 1
                continue
            if r.is_duplicate:
                duplicate += 1
                continue
            if r.is_secondary:
                secondary += 1
                continue
            if r.is_unmapped:
                unmapped += 1
                continue
            if not _is_read_uniq_mapping(r):
                multi += 1
                continue

            map_strand = '-' if r.is_reverse else '+'
            ref_positions = r.get_reference_positions()
            strand = None
            pos = None
            chrom = r.reference_name
            length = r.query_length
            if protocol == 'forward':
                if map_strand == '+':
                    strand = '+'
                    pos = ref_positions[0]
                else:
                    strand = '-'
                    pos = ref_positions[-1]
            elif protocol == 'reverse':
                if map_strand == '+':
                    strand = '-'
                    pos = ref_positions[-1]
                else:
                    strand = '+'
                    pos = ref_positions[0]
            alignments[length][strand][(chrom, pos)] += 1
            read_lengths[length] += 1

            valid += 1

            pbar.update()

    summary = ('summary:\n\ttotal_reads: {}\n\tunique_mapped: {}\n'
               '\tqcfail: {}\n\tduplicate: {}\n\tsecondary: {}\n'
               '\tunmapped:{}\n\tmulti:{}\n\nlength dist:\n').format(
                   total_count, valid, qcfail, duplicate, secondary, unmapped,
                   multi)

    for length in read_lengths:
        summary += '\t{}: {}\n'.format(length, read_lengths[length])

    with open('{}_bam_summary.txt'.format(prefix), 'w') as output:
        output.write(summary)

    return (alignments, read_lengths)


def align_metagenes(metagenes, read_lengths, prefix):
    """align metagene coverages to determine the lag of the psites

    Parameters
    ----------
    metagenes: dict
               key is the length, value is the metagene coverage
    read_lengths: dict
                  key is the length, value is the number of reads
    prefix: str
            prefix for output files

    Returns
    -------
    psite_offsets: dict
                   key is the length, value is the offset
    """
    base = int(base)
    with open(coverages) as f:
        cov_lens = f.readlines()
    cov_lens = {
        int(x.strip().split()[0]): x.strip().split()[1]
        for x in cov_lens
    }
    if base not in cov_lens:
        print('Failed to find base {} in coverages.'.format(base))
        return
    reference = pd.read_table(cov_lens[base])['count']
    to_write = 'relative lag to base: {}\n'.format(base)
    for length, path in cov_lens.items():
        cov = pd.read_table(path)['count']
        xcorr = np.correlate(reference, cov, 'full')
        origin = len(xcorr) // 2
        bound = min(base, length)
        xcorr = xcorr[origin - bound:origin + bound]
        lag = np.argmax(xcorr) - len(xcorr) // 2
        to_write += '\tlag of {}: {}\n'.format(length, lag)
    with open(saveto, 'w') as output:
        output.write(to_write)


def merge_lengths(alignments, read_lengths, psite_offsets):
    """
    Parameters
    ----------
    alignments: dict(dict(Counter))
                bam split by length, strand
    read_lengths: dict
                  key is the length, value is the number of reads
    psite_offsets: dict
                   key is the length, value is the offset
    Returns
    -------
    merged_alignments: dict(dict)
                       alignments by merging all lengths
    """
    pass


def parse_annotation(annotation):
    """
    Parameters
    ----------
    annotation: string
          path of annotation file generated by prepare_orfs
    
    Returns
    -------
    cds: List[PutativeORF]
         list of cds
    uorfs: List[PutativeORF]
          list of putative ORFs from 5'UTR
    dorfs: List[PutativeORF]
          list of putative ORFs from 3'UTR
    """

    cds = []
    uorfs = []
    dorfs = []
    print('parsing putative ORFs...')
    with open(annotation, 'r') as anno:
        total_lines = len(['' for line in anno])
    with open(annotation, 'r') as anno:
        with tqdm(total=total_lines) as pbar:
            header = True
            for line in anno:
                if header:
                    header = False
                    continue
                orf = PutativeORF.from_string(line)
                if orf is None:
                    continue
                if orf.category == 'CDS':
                    cds.append(orf)
                elif orf.category == 'uORF':
                    uorfs.append(orf)
                elif orf.category == 'dORF':
                    dorfs.append(orf)
                pbar.update()
    return (cds, uorfs, dorfs)


def orf_coverage_length(orf, alignments, length, offset_5p=0, offset_3p=0):
    """
    Parameters
    ----------
    orf: PutativeORF
         instance of PutativeORF
    alignments: dict(dict(Counter))
                alignments summarized from bam
    length: int
            the target length
    offset_5p: int
               the number of nts to include from 5'prime
    offset_3p: int
               the number of nts to include from 3'prime

    Returns
    -------
    coverage: Series
              coverage for ORF for specific length
    """
    coverage = []
    chrom = orf.chrom
    strand = orf.strand
    if strand == '-':
        offset_5p, offset_3p = offset_3p, offset_5p
    first, last = orf.intervals[0], orf.intervals[-1]
    for pos in range(first.start - offset_5p, first.start):
        try:
            coverage.append(alignments[length][strand][(chrom, pos)])
        except KeyError:
            coverage.append(0)

    for iv in orf.intervals:
        for pos in range(iv.start, iv.end + 1):
            try:
                coverage.append(alignments[length][strand][(chrom, pos)])
            except KeyError:
                coverage.append(0)

    for pos in range(last.end + 1, last.end + offset_3p + 1):
        try:
            coverage.append(alignments[length][strand][(chrom, pos)])
        except KeyError:
            coverage.append(0)

    if strand == '-':
        coverage.reverse()
        return pd.Series(
            np.array(coverage),
            index=np.arange(-offset_3p,
                            len(coverage) - offset_3p))
    else:
        return pd.Series(
            np.array(coverage),
            index=np.arange(-offset_5p,
                            len(coverage) - offset_5p))


def metagene_coverage(cds,
                      alignments,
                      prefix,
                      max_positions=500,
                      offset_5p=0,
                      offset_3p=0,
                      alignby='start_codon'):
    """
    Parameters
    ----------
    cds: List[PutativeORF]
         list of cds
    alignments: dict(dict(Counter))
                alignments summarized from bam
    prefix: str
            prefix for the output file
    max_positions: int
                   the number of nts to include
    offset_5p: int
               the number of nts to include from the 5'prime
    offset_3p: int
               the number of nts to include from the 3'prime
    alignby: str
             'start_codon' or 'stop_codon'
             align gene coverage at start or stop codon to generate metagene
             coverage

    Returns
    -------
    metagenes: dict
               key is the length, value is the metagene coverage
    """
    metagenes = {}
    for length in alignments.keys():

        metagene_coverage = pd.Series()

        for orf in cds:
            coverage = orf_coverage_length(orf, alignments, length, offset_5p,
                                           offset_3p)
            if len(coverage.index) > 0:
                min_index = min(coverage.index.tolist())
                max_index = max(coverage.index.tolist())
                coverage = coverage[np.arange(min_index,
                                              min(max_index, max_positions))]
            if coverage.mean() > 0:
                metagene_coverage = metagene_coverage.add(
                    coverage, fill_value=0)
        metagenes[length] = metagene_coverage

    return metagenes


def plot_read_lengths(read_lengths, prefix):
    """
    Parameters
    ----------
    read_lengths: dict
                  key is the length, value is the number of reads
    prefix: str
            prefix for the output file
    """
    pass


def plot_metagene(metagenes, read_lengths, prefix, offset=60):
    """
    Parameters
    ----------
    metagenes: dict
               key is the length, value is the metagene coverage
    read_lengths: dict
                  key is the length, value is the number of reads
    prefix: str
            prefix for the output file
    """
    total_reads = sum(read_lengths.values())
    with PdfPages('{}_metagene_plots.pdf'.format(prefix)) as pdf:
        for length in metagenes:
            metagene_cov = metagenes[length]
            min_index = min(metagene_cov.index.tolist())
            max_index = max(metagene_cov.index.tolist())
            offset = min(offset, max_index)
            metagene_cov = metagene_cov[np.arange(min_index, offset)]
            x = np.arange(min_index, offset)
            colors = np.tile(['r', 'g', 'b'], len(x)//3 + 1)
            xticks = np.arange(min_index, offset, 20)
            ratio = '{:.2%}'.format(read_lengths[length] / total_reads)
            fig, ax = plt.subplots()
            ax.vlines(x, ymin=np.zeros(len(x)), ymax=metagene_cov)
            ax.tick_params(axis='x', which='both', top='off', direction='out')
            ax.set_xticks(xticks)
            ax.set_xlim((min_index, offset))
            ax.set_xlabel('Distance from start codon (nt)')
            ax.set_ylabel('Number of reads')
            ax.set_title('{} nt reads, proportion: {}'.format(length, ratio))

            fig.tight_layout()
            pdf.savefig(fig)
            plt.close()


def export_orf_coverages(orfs, merged_alignments, prefix):
    """
    Parameters
    ----------
    orfs: List[PutativeORF]
          a list of putative orfs
    merged_alignments: dict(dict)
                       alignments by merging all lengths
    prefix: str
            prefix for output file
    """
    pass


def export_wig(merged_alignments, prefix):
    """
    Parameters
    ----------
    merged_alignments: dict(dict)
                       alignments by merging all lengths
    prefix: str
            prefix of output wig files
    """
    pass


def detect_orfs(gtf, fasta, bam, prefix, annotation=None, protocol=None):
    """
    Parameters
    ----------
    gtf: str
         Path to the GTF file
    fasta: str
           Path to the FASTA file
    bam: str
         Path to the bam file
    prefix: str
            prefix for all output files
    annotation: str
                Path for annontation files of putative ORFs
                It will be automatically generated if None
    protocol: str
              'forward' for stranded, 'reverse' for reverse stranded
              It will be automatically inferred if None
    """

    if not isinstance(gtf, GTFReader):
        gtf = GTFReader(gtf)

    if not isinstance(fasta, FastaReader):
        fasta = FastaReader(fasta)

    if annotation is None:
        cds, uorfs, dorfs = prepare_orfs(gtf, fasta, prefix)
    else:
        cds, uorfs, dorfs = parse_annotation(annotation)

    if protocol is None:
        protocol = infer_protocol(bam, gtf, prefix)

    alignments, read_lengths = split_bam(bam, protocol, prefix)
    plot_read_lengths(read_lengths, prefix)
    metagenes = metagene_coverage(cds, alignments, prefix)
    plot_metagene(metagenes, read_lengths, prefix)
    psite_offsets = align_metagenes(metagenes, read_lengths, prefix)
    merged_alignments = merge_lengths(alignments, read_lengths, psite_offsets)
    export_wig(merged_alignments, prefix)
    export_orf_coverages(cds + uorfs + dorfs, merged_alignments, prefix)
