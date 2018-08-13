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
from .interval import Interval
from .count import _is_read_uniq_mapping

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
