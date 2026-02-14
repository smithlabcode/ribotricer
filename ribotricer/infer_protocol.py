"""Infer experimental protocol"""

# Part of ribotricer software
#
# Copyright (C) 2020 Saket Choudhary, Wenzheng Li, and Andrew D Smith
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

from __future__ import annotations

from collections import Counter, defaultdict
from typing import TYPE_CHECKING

import pysam
from quicksect import Interval, IntervalTree

from .common import is_read_uniq_mapping

if TYPE_CHECKING:
    pass

# required to convert numeric strands to '+/-'
NUM_TO_STRAND: dict[int, str] = {1: "+", -1: "-"}


def infer_protocol(
    bam: str,
    gene_interval_tree: defaultdict[str, IntervalTree],
    prefix: str,
    n_reads: int = 20000,
) -> str:
    """Infer strandedness protocol given a bam file.

    Parameters
    ----------
    bam : str
        Path to bam file.
    gene_interval_tree : defaultdict[str, IntervalTree]
        Dictionary mapping chrom to (start, end, strand) interval tree.
    prefix : str
        Prefix for protocol file.
    n_reads : int, optional
        Number of reads to use (downsampled), by default 20000.

    Returns
    -------
    str
        Protocol string: 'forward' or 'reverse'.

    Notes
    -----
    The strategy to do this is simple: keep a track
    of mapped reads and their strand and then tally
    if the location of their mapping has a gene defined
    on the positive strand or the negative strand.

    If the first and second characters denote the mapping and
    gene strand respectively:
    Higher proportion of (++, --) implies forward protocol
    Higher proportion of (+-, -+) implies reverse protocol
    Equal proportion of the above two scenarios implies unstranded protocol.
    """
    iteration = 0
    bam_file = pysam.AlignmentFile(bam, "rb")
    strandedness: Counter[str] = Counter()

    for read in bam_file.fetch(until_eof=True):
        if iteration <= n_reads:
            if is_read_uniq_mapping(read):
                if read.is_reverse:
                    mapped_strand = "-"
                else:
                    mapped_strand = "+"
                mapped_start = read.reference_start
                mapped_end = read.reference_end
                chrom = read.reference_name

                if chrom is not None and mapped_end is not None:
                    # get corresponding gene's strand
                    interval = list(
                        set(
                            gene_interval_tree[chrom].find(
                                Interval(mapped_start, mapped_end)
                            )
                        )
                    )
                    if len(interval) == 1:
                        # Filter out genes with ambiguous strand info
                        # (those) that have a tx_start on opposite strands
                        gene_strand = NUM_TO_STRAND[interval[0].data]
                        # count table for mapped strand vs gene strand
                        strandedness["{}{}".format(mapped_strand, gene_strand)] += 1
                        iteration += 1

    bam_file.close()

    # Add pseudocounts
    strandedness["++"] += 1
    strandedness["--"] += 1
    strandedness["+-"] += 1
    strandedness["-+"] += 1

    total = sum(strandedness.values())
    forward_mapped_reads = strandedness["++"] + strandedness["--"]
    reverse_mapped_reads = strandedness["-+"] + strandedness["+-"]
    to_write = (
        "In total {} reads checked:\n"
        '\tNumber of reads explained by "++, --": {} ({:.4f})\n'
        '\tNumber of reads explained by "+-, -+": {} ({:.4f})\n'
    ).format(
        total,
        forward_mapped_reads,
        forward_mapped_reads / total,
        reverse_mapped_reads,
        reverse_mapped_reads / total,
    )
    with open("{}_protocol.txt".format(prefix), "w") as output:
        output.write(to_write)
    protocol = "forward"
    if reverse_mapped_reads > forward_mapped_reads:
        protocol = "reverse"
    return protocol
