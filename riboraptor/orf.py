"""Utilities for translating ORF detection
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import warnings

from collections import Counter
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam

from .wig import WigReader
from .fasta import FastaReader
from .interval import Interval
from .count import _is_read_uniq_mapping
from .helpers import merge_intervals

def transcript_to_genome_iv(start, end, intervals, reverse=False):
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


def search_orfs(fasta, intervals):
    if not intervals:
        return []

    orfs = []
    if not isinstance(fasta, FastaReader):
        fasta = FastaReader(fasta)
    intervals, _, _ = merge_intervals(intervals)
    sequences = fasta.query(intervals)
    merged_seq = ''.join(sequences)
    reverse = False
    strand = intervals[0].strand
    if strand == '-':
        merged_seq = merged_seq[::-1]
        reverse = True
    start_codons = ['ATG', 'CTG', 'GTG']
    stop_codons = ['TAG', 'TAA', 'TGA']
    for sc in start_codons:
        cur = 0
        while cur < len(merged_seq):
            start = merged_seq.find(sc, cur)
            if start == -1:
                break
            cur = start + 1
            for i in range(start, len(merged_seq), 3):
                if merged_seq[i:i+3] in stop_codons:
                    ### found orf
                    ivs = transcript_to_genome_iv(start, i+2, intervals, reverse)
                    if ivs:
                        orfs.append(ivs)
                    break
    return orfs


def prepare_orfs(gtf, fasta, prefix):

    gtf = pd.read_table(gtf, header=None, skiprows=5)
    gtf.rename(columns={0: 'seqname', 1: 'source', 2: 'feature', 3: 'start',
                        4: 'end', 5: 'score', 6: 'strand', 7: 'frame',
                        8: 'attribute'}, inplace=True)
    gtf = gtf[gtf['feature'].isin(['CDS', 'UTR'])]
    cds = gtf[gtf['feature'] == 'CDS']
    utr = gtf[gtf['feature'] == 'UTR']

    iteration = 0
    ### process CDS gtf
    cds_intervals = defaultdict(lambda: defaultdict(list))
    for i, r in cds.iterrows():
        iteration += 1
        if iteration % 1000 == 0:
            print('{} gtf records processed.'.format(iteration))
        chrom = r['seqname']
        feature = r['feature']
        start = r['start']
        end = r['end']
        strand = r['strand']
        attribute = r['attribute']
        attribute = {x.strip().split()[0]: x.strip().split()[1]
                        for x in attribute.split(';') if len(x.split()) > 1}
        if 'transcript_id' not in attribute:
            print('missing transcript_id {}: {}-{}'.format(chrom, start, end))
            continue
        transcript_id = attribute['transcript_id']
        if 'gene_id' not in attribute:
            print('missing gene_id {}: {}-{}'.format(chrom, start, end))
            continue
        gene_id = attribute['gene_id']
        cds_intervals[gene_id][transcript_id].append(Interval(chrom, start, end,
            strand))
    
    ### process UTR gtf
    utr5_intervals = defaultdict(list)
    utr3_intervals = defaultdict(list)
    for i, r in utr.iterrows():
        iteration += 1
        if iteration % 1000 == 0:
            print('{} gtf records processed.'.format(iteration))
        chrom = r['seqname']
        start = r['start']
        end = r['end']
        strand = r['strand']
        attribute = r['attribute']
        attribute = {x.strip().split()[0]: x.strip().split()[1]
                        for x in attribute.split(';') if len(x.split()) > 1}
        if 'transcript_id' not in attribute:
            print('missing transcript_id {}: {}-{}'.format(chrom, start, end))
            continue
        transcript_id = attribute['transcript_id']
        if 'gene_id' not in attribute:
            print('missing gene_id {}: {}-{}'.format(chrom, start, end))
            continue
        gene_id = attribute['gene_id']
        if (gene_id not in cds_intervals or transcript_id not in
                cds_intervals[gene_id]):
            print('missing CDS for UTR {}: {}-{}'.format(chrom, start, end))
            continue
        gene_cds = []
        for transcript in cds_intervals[gene_id]:
            gene_cds += cds_intervals[gene_id][transcript]
        first_cds = gene_cds[0]
        for gc in gene_cds:
            if gc.start < first_cds.start:
                first_cds = gc
        last_cds = gene_cds[-1]
        for gc in gene_cds:
            if gc.end > last_cds.end:
                last_cds = gc
        interval = Interval(chrom, start, end, strand)
        if start < first_cds.start:
            if end >= first_cds.start:
                interval.end = first_cds.start - 1
            if strand == '+':
                utr5_intervals[transcript_id].append(interval)
            else:
                utr3_intervals[transcript_id].append(interval)
        elif end > last_cds.end:
            if start <= last_cds.end:
                interval.start = last_cds.end + 1
            if strand == '+':
                utr3_intervals[transcript_id].append(interval)
            else:
                utr5_intervals[transcript_id].append(interval)

    ### sort the intervals
    for transcript in utr5_intervals:
        utr5_intervals[transcript].sort(key=lambda x: x.start)
    for transcript in utr3_intervals:
        utr3_intervals[transcript].sort(key=lambda x: x.start)

    if not isinstance(fasta, FastaReader):
        fasta = FastaReader(fasta)

    uorfs = {}
    print('searching uORFs...')
    for tid, invs in utr5_intervals.items():
        orfs = search_orfs(fasta, invs)
        for orf in orfs:
            start = orf[0].start
            end = orf[-1].end
            total_len = sum(i.end - i.start + 1 for i in orf)
            orf_id = '{}_{}_{}_{}'.format(tid, start, end, total_len)
            uorfs[orf_id] = orf

    dorfs = {}
    print('searching dORFs...')
    for tid, invs in utr3_intervals.items():
        orfs = search_orfs(fasta, invs)
        for orf in orfs:
            start = orf[0].start
            end = orf[-1].end
            total_len = sum(i.end - i.start + 1 for i in orf)
            orf_id = '{}_{}_{}_{}'.format(tid, start, end, total_len)
            dorfs[orf_id] = orf

    ### save to bed file
    cds_bed = ''
    for gid in cds_intervals:
        for tid in cds_intervals[gid]:
            for iv in cds_intervals[gid][tid]:
                cds_bed += '{}\t{}\t{}\t{}\t{}\t{}\n'.format(iv.chrom, iv.start,
                        iv.end, tid, '.', iv.strand)
    with open('{}_cds.bed'.format(prefix), 'w') as output:
        output.write(cds_bed)

    utr5_bed = ''
    for oid in uorfs:
        for iv in uorfs[oid]:
            utr5_bed += '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(iv.chrom,
                    iv.start, iv.end, oid, '.', iv.strand)
    with open('{}_utr5.bed'.format(prefix), 'w') as output:
        output.write(utr5_bed)

    utr3_bed = ''
    for oid in dorfs:
        for iv in dorfs[oid]:
            utr3_bed += '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(iv.chrom,
                    iv.start, iv.end, oid, '.', iv.strand)
    with open('{}_utr3.bed'.format(prefix), 'w') as output:
        output.write(utr3_bed)


def split_bam(bam, protocol, prefix):
    """Split bam by read length and strand

    Parameters
    ----------
    bam : str
          Path to bam file
    protocol: str
          Experiment protocol [forward, reverse]
    prefix: str
            prefix for output files: {prefix}_xxnt_pos.wig and
            {prefix}__xxnt_neg.wig
    """
    coverages = defaultdict(lambda: defaultdict(Counter))
    iteration=qcfail=duplicate=secondary=unmapped=multi=valid=0
    bam = pysam.AlignmentFile(bam, 'rb')
    for r in bam.fetch(until_eof=True):

        iteration += 1
        if iteration % 1000 == 0:
            print('{} reads processed.'.format(iteration))

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
                pos =  ref_positions[-1]
            else:
                strand = '+'
                pos = ref_positions[0]
        coverages[length][strand][(chrom, pos)] += 1

        valid += 1
        
    summary = 'summary:\n\ttotal_reads: {}\n\tunique_mapped: {}\n' \
              '\tqcfail: {}\n\tduplicate: {}\n\tsecondary: {}\n' \
              '\tunmapped:{}\n\tmulti:{}\n\nlength dist:\n'.format(iteration,
                       valid, qcfail, duplicate, secondary, unmapped, multi)

    for length in coverages:
        reads_of_length = 0
        for strand in coverages[length]:
            to_write = ''
            cur_chrom = ''
            for chrom, pos in sorted(coverages[length][strand]):
                if chrom != cur_chrom:
                    cur_chrom = chrom
                    to_write += 'variableStep chrom={}\n'.format(chrom)
                to_write += '{}\t{}\n'.format(pos,
                        coverages[length][strand][(chrom, pos)])
            fname = '{}_{}nt_{}.wig'.format(prefix, length,
                    'pos' if strand == '+' else 'neg')
            with open(fname, 'w') as output:
                output.write(to_write)
            reads_of_length += 1
        summary += '\t{}: {}\n'.format(length, reads_of_length)
    with open('{}_summary.txt'.format(prefix), 'w') as output:
        output.write(summary)


def align_coverages(coverages, base, saveto):
    """align coverages to determine the lag to the base

    Parameters
    ----------
    coverages: str
               Path to file which contains paths of all metagene
               from different lengths
               format:
               length (e.g. 28) path (e.g. metagene_28.tsv)
               length (e.g. 29) path (e.g. metagene_29.tsv)
    base: int
          The reference length to align against
    saveto: str
          Path to save the aligned offsets
    """
    base = int(base)
    with open(coverages) as f:
        cov_lens = f.readlines()
    cov_lens = {int(x.strip().split()[0]): x.strip().split()[1]
                    for x in cov_lens}
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
        xcorr = xcorr[origin-bound: origin+bound]
        lag = np.argmax(xcorr) - len(xcorr)//2
        to_write += '\tlag of {}: {}\n'.format(length, lag)
    with open(saveto, 'w') as output:
        output.write(to_write)

def merge_wigs(wigs, offsets, strand, saveto):
    """merge wigs from different lengths into one with shift of offsets

    Parameters
    ----------
    wigs: str
          Path to file which contains paths of all wigs from differnt lengths
          format:
          length1 path1
          length2 path2
    offsets: str
             Path to file which contains offset for each length
             format:
             length1 offset1
             length2 offset2
    strand: str
            '+' for positive strand,
            '-' for negative strand
    saveto: str
            Path to save merged wig
    """
    coverages = defaultdict(int)
    with open(wigs) as wf:
        wigs = {int(x.strip().split()[0]): x.strip().split()[1]
                    for x in wf.readlines()}
    with open(offsets) as of:
        offsets = {int(x.strip().split()[0]): int(x.strip().split()[1])
                       for x in of.readlines()}
    for length, wig in wigs.items():
        with open(wig) as f:
            for line in f:
                if line.startswith('variableStep'):
                    line = line.strip()
                    chrom = line[line.index('=')+1:]
                else:
                    pos, count = line.strip().split()
                    pos, count = int(pos), int(count)
                    if strand == '+':
                        pos_shifted = pos + offsets[length]
                    else:
                        pos_shifted = pos - offsets[length]
                    if pos_shifted >= 0:
                        coverages[(chrom, pos_shifted)] += count
    to_write = ''
    cur_chrom = ''
    for chrom, pos in sorted(coverages):
        if chrom != cur_chrom:
            cur_chrom = chrom
            to_write += 'variableStep chrom={}\n'.format(chrom)
        to_write += '{}\t{}\n'.format(pos, coverages[(chrom, pos)])
    with open(saveto, 'w') as output:
        output.write(to_write)
