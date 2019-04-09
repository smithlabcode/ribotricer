"""Utilities for spliting bam file"""
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

import warnings

from collections import Counter
from collections import defaultdict

import pysam
from tqdm import *

from .const import TYPICAL_OFFSET
from .common import is_read_uniq_mapping


def split_bam(bam, protocol, prefix, read_lengths=None):
    """Split bam by read length and strand

    Parameters
    ----------
    bam : str
          Path to bam file
    protocol: str
          Experiment protocol [forward, reverse]
    prefix: str
            prefix for output files
    read_lengths: list[int]
                  read lengths to use
                  If None, it will be automatically determined by assessing
                  the periodicity of metagene profile of this read length

    Returns
    -------
    alignments: dict(dict(Counter))
                bam split by length, strand, (chrom, pos)
    read_length_counts: dict
                  key is the length, value is the number of reads
    """
    alignments = defaultdict(lambda: defaultdict(Counter))
    read_length_counts = defaultdict(int)
    total_count = qcfail = duplicate = secondary = unmapped = multi = valid = 0
    # print('reading bam file...')
    bam = pysam.AlignmentFile(bam, 'rb')
    for r in tqdm(bam.fetch(until_eof=True)):

        total_count += 1

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
        if not is_read_uniq_mapping(r):
            multi += 1
            continue

        map_strand = '-' if r.is_reverse else '+'
        ref_positions = r.get_reference_positions()
        strand = None
        pos = None
        chrom = r.reference_name
        # length = r.query_length
        length = len(ref_positions)
        if read_lengths is not None and length not in read_lengths:
            continue
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

        # convert bam coordinate to one-based
        alignments[length][strand][(chrom, pos + 1)] += 1
        read_length_counts[length] += 1

        valid += 1

    summary = ('summary:\n\ttotal_reads: {}\n\tunique_mapped: {}\n'
               '\tqcfail: {}\n\tduplicate: {}\n\tsecondary: {}\n'
               '\tunmapped:{}\n\tmulti:{}\n\nlength dist:\n').format(
                   total_count, valid, qcfail, duplicate, secondary, unmapped,
                   multi)

    for length in sorted(read_length_counts):
        summary += '\t{}: {}\n'.format(length, read_length_counts[length])

    with open('{}_bam_summary.txt'.format(prefix), 'w') as output:
        output.write(summary)

    return (alignments, read_length_counts)
