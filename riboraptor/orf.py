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
    bam = pysam.AlignmentFile(bam)
    for r in bam:
        if r.is_qcfail: continue
        if r.is_duplicate: continue
        if r.is_secondary: continue
        if r.is_unmapped: continue
        if not _is_read_uniq_mapping(r): continue

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

    for length in coverages:
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
