"""Utilities for spliting bam file"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import warnings

from collections import Counter
from collections import defaultdict

import pysam
from tqdm import *

from .common import is_read_uniq_mapping


def split_bam(bam, protocol, prefix):
    """Split bam by read length and strand

    Parameters
    ----------
    bam : str
          Path to bam file
    protocol: str
          Experiment protocol [forward, reverse]
    prefix: str
            prefix for output files

    Returns
    -------
    alignments: dict(dict(Counter))
                bam split by length, strand, (chrom, pos)
    read_lengths: dict
                  key is the length, value is the number of reads
    """
    alignments = defaultdict(lambda: defaultdict(Counter))
    read_lengths = defaultdict(int)
    total_count = qcfail = duplicate = secondary = unmapped = multi = valid = 0
    print('reading bam file...')
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
        read_lengths[length] += 1

        valid += 1

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
