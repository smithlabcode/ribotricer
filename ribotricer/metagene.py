"""Metagene profile related functions"""
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

import sys
from collections import Counter, OrderedDict

import numpy as np
import pandas as pd
from tqdm import tqdm

from .const import CUTOFF, TYPICAL_OFFSET
from .interval import Interval
from .statistics import coherence


def next_genome_pos(ivs, max_positions, leader, trailer, reverse=False):
    if len(ivs) == 0:
        return iter([])
    cnt = 0
    leader_iv = Interval(
        ivs[0].chrom, ivs[0].start - leader, ivs[0].start - 1, ivs[0].strand
    )
    trailer_iv = Interval(
        ivs[-1].chrom, ivs[-1].end + 1, ivs[-1].end + trailer, ivs[-1].strand
    )
    combined_ivs = [leader_iv] + ivs + [trailer_iv]

    cnt = 0
    if not reverse:
        for interval in combined_ivs:
            for pos in range(interval.start, interval.end + 1):
                cnt += 1
                if cnt <= max_positions:
                    yield pos
    else:
        for interval in reversed(combined_ivs):
            for pos in range(interval.end, interval.start - 1, -1):
                cnt += 1
                if cnt <= max_positions:
                    yield pos


def orf_coverage_length(
    orf, alignments, length, max_positions, offset_5p=20, offset_3p=0
):
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
    if strand == "-":
        offset_5p, offset_3p = offset_3p, offset_5p

    for pos in next_genome_pos(
        orf.intervals, max_positions, offset_5p, offset_3p, strand == "-"
    ):
        try:
            coverage.append(alignments[length][strand][(chrom, pos)])
        except KeyError:
            coverage.append(0)

    if strand == "-":
        from_start = pd.Series(
            np.array(coverage), index=np.arange(-offset_3p, len(coverage) - offset_3p)
        )
        from_stop = pd.Series(
            np.array(coverage),
            index=np.arange(offset_5p - len(coverage) + 1, offset_5p + 1),
        )
    else:
        from_start = pd.Series(
            np.array(coverage), index=np.arange(-offset_5p, len(coverage) - offset_5p)
        )
        from_stop = pd.Series(
            np.array(coverage),
            index=np.arange(offset_3p - len(coverage) + 1, offset_3p + 1),
        )

    return (from_start, from_stop)


def metagene_coverage(
    cds,
    alignments,
    read_lengths,
    prefix,
    max_positions=600,
    offset_5p=20,
    offset_3p=0,
    meta_min_reads=100000,
):
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
    # print('calculating metagene profiles...')
    metagenes = {}

    # remove read length whose read number is small
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
                orf, alignments, length, max_positions, offset_5p, offset_3p
            )
            cov_mean = from_start.mean()
            if cov_mean > 0:
                from_start = from_start / cov_mean
                from_start = from_start.fillna(0)
                metagene_coverage_start = metagene_coverage_start.add(
                    from_start, fill_value=0
                )
                position_counter_start += Counter(from_start.index.tolist())

                from_stop = from_stop / cov_mean
                from_stop = from_stop.fillna(0)
                metagene_coverage_stop = metagene_coverage_stop.add(
                    from_stop, fill_value=0
                )
                position_counter_stop += Counter(from_stop.index.tolist())

        if len(position_counter_start) != len(metagene_coverage_start) or len(
            position_counter_stop
        ) != len(metagene_coverage_stop):
            raise RuntimeError("Metagene coverage and counter mismatch")
        position_counter_start = pd.Series(position_counter_start)
        metagene_coverage_start = metagene_coverage_start.div(position_counter_start)
        position_counter_stop = pd.Series(position_counter_stop)
        metagene_coverage_stop = metagene_coverage_stop.div(position_counter_stop)

        coh_5p, valid_5p = coherence(metagene_coverage_start.tolist())
        coh_3p, valid_3p = coherence(metagene_coverage_stop.tolist())
        metagenes[length] = (
            metagene_coverage_start,
            metagene_coverage_stop,
            coh_5p,
            valid_5p,
            coh_3p,
            valid_3p,
        )

    to_write_5p = "fragment_length\toffset_5p\tprofile\tphase_score\tvalid_codons\n"
    to_write_3p = "fragment_length\toffset_3p\tprofile\tphase_score\tvalid_codons\n"
    for length in sorted(metagenes):
        to_write_5p += "{}\t{}\t{}\t{}\t{}\n".format(
            length,
            offset_5p,
            metagenes[length][0].tolist(),
            metagenes[length][2],
            metagenes[length][3],
        )
        to_write_3p += "{}\t{}\t{}\t{}\t{}\n".format(
            length,
            offset_3p,
            metagenes[length][1].tolist(),
            metagenes[length][4],
            metagenes[length][5],
        )

    with open("{}_metagene_profiles_5p.tsv".format(prefix), "w") as output:
        output.write(to_write_5p)
    with open("{}_metagene_profiles_3p.tsv".format(prefix), "w") as output:
        output.write(to_write_3p)

    return metagenes


def align_metagenes(
    metagenes, read_lengths, prefix, phase_score_cutoff=CUTOFF, remove_nonperiodic=False
):
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
    remove_nonperiodic: bool
                        Whether remove non-periodic read lengths

    Returns
    -------
    psite_offsets: dict
                   key is the length, value is the offset
    """
    # print('aligning metagene profiles from different lengths...')

    # discard non-periodic read lengths
    if remove_nonperiodic:
        for length, (_, _, coh, _, _, _) in list(metagenes.items()):
            if coh < phase_score_cutoff:
                del read_lengths[length]
                del metagenes[length]

    if len(read_lengths) == 0:
        sys.exit(
            "Warning: no periodic read length found... using cutoff {}".format(
                phase_score_cutoff
            )
        )

    psite_offsets = OrderedDict()
    base = n_reads = 0
    for length, reads in list(read_lengths.items()):
        if reads > n_reads:
            base = length
            n_reads = reads
    reference = metagenes[base][0].values
    to_write = "relative lag to base: {}\n".format(base)
    for length, (meta, _, _, _, _, _) in list(metagenes.items()):
        cov = meta.values
        xcorr = np.correlate(reference, cov, "full")
        origin = len(xcorr) // 2
        bound = min(base, length)
        xcorr = xcorr[origin - bound : origin + bound]
        lag = np.argmax(xcorr) - len(xcorr) // 2
        psite_offsets[length] = lag + TYPICAL_OFFSET
        to_write += "\tlag of {}: {}\n".format(length, lag)
    with open("{}_psite_offsets.txt".format(prefix), "w") as output:
        output.write(to_write)
    return psite_offsets
