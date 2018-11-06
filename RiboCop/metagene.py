"""Metagene profile related functions"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import warnings

from collections import Counter
from collections import defaultdict

import sys
from tqdm import *
import numpy as np
import pandas as pd

from .interval import Interval
from .statistics import coherence
from .constants import CUTOFF


def next_genome_pos(ivs, max_positions, leader, trailer, reverse=False):
    if len(ivs) == 0:
        return iter([])
    cnt = 0
    leader_iv = Interval(ivs[0].chrom, ivs[0].start - leader, ivs[0].start - 1,
                         ivs[0].strand)
    trailer_iv = Interval(ivs[-1].chrom, ivs[-1].end + 1,
                          ivs[-1].end + trailer, ivs[-1].strand)
    combined_ivs = [leader_iv] + ivs + [trailer_iv]

    cnt = 0
    if not reverse:
        for iv in combined_ivs:
            for pos in range(iv.start, iv.end + 1):
                cnt += 1
                if cnt > max_positions:
                    break
                yield pos
    else:
        for iv in reversed(combined_ivs):
            for pos in range(iv.end, iv.start - 1, -1):
                cnt += 1
                if cnt > max_positions:
                    break
                yield pos


def orf_coverage_length(orf,
                        alignments,
                        length,
                        max_positions,
                        offset_5p=20,
                        offset_3p=0):
    """
    Parameters
    ----------
    orf: ORF
         instance of ORF
    alignments: dict(dict(Counter))
                alignments summarized from bam
    length: int
            the target length
    max_positions: int
                   the number of nts to include
    offset_5p: int
               the number of nts to include from 5'prime
    offset_3p: int
               the number of nts to include from 3'prime

    Returns
    -------
    from_start: Series
                coverage for ORF for specific length aligned at start codon
    from_stop: Series
               coverage for ORF for specific length aligned at stop codon
    """
    coverage = []
    chrom = orf.chrom
    strand = orf.strand
    if strand == '-':
        offset_5p, offset_3p = offset_3p, offset_5p

    for pos in next_genome_pos(orf.intervals, max_positions, offset_5p,
                               offset_3p, strand == '-'):
        try:
            coverage.append(alignments[length][strand][(chrom, pos)])
        except KeyError:
            coverage.append(0)

    if strand == '-':
        from_start = pd.Series(
            np.array(coverage),
            index=np.arange(-offset_3p,
                            len(coverage) - offset_3p))
        from_stop = pd.Series(
            np.array(coverage),
            index=np.arange(offset_5p - len(coverage) + 1, offset_5p + 1))
    else:
        from_start = pd.Series(
            np.array(coverage),
            index=np.arange(-offset_5p,
                            len(coverage) - offset_5p))
        from_stop = pd.Series(
            np.array(coverage),
            index=np.arange(offset_3p - len(coverage) + 1, offset_3p + 1))

    return (from_start, from_stop)


def metagene_coverage(cds,
                      alignments,
                      read_lengths,
                      prefix,
                      max_positions=600,
                      offset_5p=20,
                      offset_3p=0,
                      meta_min_reads=100000):
    """
    Parameters
    ----------
    cds: List[ORF]
         list of cds
    alignments: dict(dict(Counter))
                alignments summarized from bam
    read_lengths: dict
                  key is the length, value is the number reads
    prefix: str
            prefix for the output file
    max_positions: int
                   the number of nts to include
    offset_5p: int
               the number of nts to include from the 5'prime
    offset_3p: int
               the number of nts to include from the 3'prime
    meta_min_reads: int
                    minimum number of reads for a read length to be considered

    Returns
    -------
    metagenes: dict
               key is the length, value is (from_start, from_stop, coherence,
               pval)
    """
    print('calculating metagene profiles...')
    metagenes = {}

    ### remove read length whose read number is small
    for length, reads in list(read_lengths.items()):
        if reads < meta_min_reads:
            del read_lengths[length]

    for length in tqdm(read_lengths):

        metagene_coverage_start = pd.Series()
        position_counter_start = Counter()
        metagene_coverage_stop = pd.Series()
        position_counter_stop = Counter()

        for orf in tqdm(cds):
            from_start, from_stop = orf_coverage_length(
                orf, alignments, length, max_positions, offset_5p, offset_3p)
            cov_mean = from_start.mean()
            if cov_mean > 0:
                from_start = from_start / cov_mean
                from_start = from_start.fillna(0)
                metagene_coverage_start = metagene_coverage_start.add(
                    from_start, fill_value=0)
                position_counter_start += Counter(from_start.index.tolist())

                from_stop = from_stop / cov_mean
                from_stop = from_stop.fillna(0)
                metagene_coverage_stop = metagene_coverage_stop.add(
                    from_stop, fill_value=0)
                position_counter_stop += Counter(from_stop.index.tolist())

        if (len(position_counter_start) != len(metagene_coverage_start)
                or len(position_counter_stop) != len(metagene_coverage_stop)):
            raise RuntimeError('Metagene coverage and counter mismatch')
        position_counter_start = pd.Series(position_counter_start)
        metagene_coverage_start = metagene_coverage_start.div(
            position_counter_start)
        position_counter_stop = pd.Series(position_counter_stop)
        metagene_coverage_stop = metagene_coverage_stop.div(
            position_counter_stop)

        coh, pval, valid = coherence(metagene_coverage_start.tolist())
        metagenes[length] = (metagene_coverage_start, metagene_coverage_stop,
                             coh, pval, valid)

    to_write = ''
    for length in sorted(metagenes):
        to_write += '{}\t{}\n'.format(
            length, metagenes[length][0].astype(int).tolist())

    with open('{}_metagene_profiles.tsv'.format(prefix), 'w') as output:
        output.write(to_write)

    return metagenes


def align_metagenes(metagenes, read_lengths, prefix):
    """align metagene coverages to determine the lag of the psites, the
    non-periodic read length will be discarded in this step

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
    print('aligning metagene profiles from different lengths...')

    ### discard non-periodic read lengths
    for length, (_, _, coh, _, _) in list(metagenes.items()):
        if coh < CUTOFF:
            del read_lengths[length]
            del metagenes[length]

    if len(read_lengths) == 0:
        print('no periodic read length found...')
        sys.exit(0)

    psite_offsets = {}
    base = n_reads = 0
    for length, reads in read_lengths.items():
        if reads > n_reads:
            base = length
            n_reads = reads
    reference = metagenes[base].values
    to_write = 'relative lag to base: {}\n'.format(base)
    for length, meta in metagenes.items():
        cov = meta.values
        xcorr = np.correlate(reference, cov, 'full')
        origin = len(xcorr) // 2
        bound = min(base, length)
        xcorr = xcorr[origin - bound:origin + bound]
        lag = np.argmax(xcorr) - len(xcorr) // 2
        psite_offsets[length] = lag
        to_write += '\tlag of {}: {}\n'.format(length, lag)
    with open('{}_psite_offsets.txt'.format(prefix), 'w') as output:
        output.write(to_write)
    return psite_offsets
