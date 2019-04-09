"""infer experimental protocol"""
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

import pysam
import numpy as np
import sys
from collections import Counter
from .gtf import GTFReader
from .common import is_read_uniq_mapping


def infer_protocol(bam, refseq, prefix, n_reads=20000):
    """Infer strandedness protocol given a bam file

    Parameters
    ----------
    bam: str
         Path to bam file
    refseq: defaultdict(IntervalTree)
            chrom: (start, end, strand)
    prefix: str
            Prefix for protocol file
    n_reads: int
             Number of reads to use (downsampled)

    Returns
    -------
    protocol: string
              forward/reverse
    """
    iteration = 0
    bam = pysam.AlignmentFile(bam, 'rb')
    strandedness = Counter()
    for read in bam.fetch(until_eof=True):
        if not is_read_uniq_mapping(read):
            continue
        if read.is_reverse:
            mapped_strand = '-'
        else:
            mapped_strand = '+'
        mapped_start = read.reference_start
        mapped_end = read.reference_end
        chrom = read.reference_name
        gene_strand = list(set(refseq[chrom].find(mapped_start, mapped_end)))
        if len(gene_strand) != 1:
            # Filter out genes with ambiguous strand info
            # (those) that have a tx_start on opposite strands
            continue
        gene_strand = gene_strand[0]
        strandedness['{}{}'.format(mapped_strand, gene_strand)] += 1
        iteration += 1
        if iteration >= n_reads:
            break
    ## Add pseudocounts
    strandedness['++'] += 1
    strandedness['--'] += 1
    strandedness['+-'] += 1
    strandedness['-+'] += 1

    total = sum(strandedness.values())
    forward_mapped_reads = (strandedness['++'] + strandedness['--'])
    reverse_mapped_reads = (strandedness['-+'] + strandedness['+-'])
    to_write = (
        'In total {} reads checked:\n'
        '\tNumber of reads explained by "++, --": {} ({:.4f})\n'
        '\tNumber of reads explained by "+-, -+": {} ({:.4f})\n').format(
            total, forward_mapped_reads, forward_mapped_reads / total,
            reverse_mapped_reads, reverse_mapped_reads / total)
    with open('{}_protocol.txt'.format(prefix), 'w') as output:
        output.write(to_write)
    protocol = 'forward'
    if reverse_mapped_reads > forward_mapped_reads:
        protocol = 'reverse'
    return protocol
