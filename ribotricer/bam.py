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

from collections import Counter
from collections import defaultdict

import pysam
from tqdm import tqdm

from .common import is_read_uniq_mapping


def split_bam(bam_path, protocol, prefix, read_lengths=None):
    """Split bam by read length and strand

    Parameters
    ----------
    bam_path : str
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
    # First pass just counts the reads
    # this is required to display a progress bar
    bam = pysam.AlignmentFile(bam_path, "rb")
    total_reads = bam.count(until_eof=True)
    bam.close()
    with tqdm(total=total_reads) as pbar:
        bam = pysam.AlignmentFile(bam_path, "rb")
        for read in bam.fetch(until_eof=True):
            pbar.update()
            # Track if the current read is usable
            is_usable = True
            total_count += 1

            if read.is_qcfail:
                qcfail += 1
                is_usable = False
            elif read.is_duplicate:
                duplicate += 1
                is_usable = False
            elif read.is_secondary:
                secondary += 1
                is_usable = False
            elif read.is_unmapped:
                unmapped += 1
                is_usable = False
            elif not is_read_uniq_mapping(read):
                multi += 1
                is_usable = False

            if is_usable:
                map_strand = "-" if read.is_reverse else "+"
                ref_positions = read.get_reference_positions()
                strand = None
                pos = None
                chrom = read.reference_name
                length = len(ref_positions)
                if read_lengths is not None and length not in read_lengths:
                    # Do nothing
                    pass
                else:
                    if protocol == "forward":
                        # Library preparation was forward-stranded:
                        # Genes defined on + strand should have
                        # reads mapping only on the positive strand
                        if map_strand == "+":
                            strand = "+"
                            # Track the 5'end
                            pos = ref_positions[0]
                        else:
                            strand = "-"
                            # For negative strand of forward protocol read the
                            # the 5'end of the read is the last element
                            pos = ref_positions[-1]
                    elif protocol == "reverse":
                        # Library preparation was reverse-stranded
                        # Mappings on the positive strand are
                        # switched to negative strand with their positions
                        # reversed and vice versa for mappings on the negative
                        # strand.
                        if map_strand == "+":
                            strand = "-"
                            # The 5' end is the last position
                            pos = ref_positions[-1]
                        else:
                            strand = "+"
                            # The 5'end is the first position
                            pos = ref_positions[0]

                    # convert bam coordinate to one-based
                    alignments[length][strand][(chrom, pos + 1)] += 1
                    read_length_counts[length] += 1
                    valid += 1
    bam.close()
    summary = (
        "summary:\n\ttotal_reads: {}\n\tunique_mapped: {}\n"
        "\tqcfail: {}\n\tduplicate: {}\n\tsecondary: {}\n"
        "\tunmapped:{}\n\tmulti:{}\n\nlength dist:\n"
    ).format(total_count, valid, qcfail, duplicate, secondary, unmapped, multi)

    for length in sorted(read_length_counts):
        summary += "\t{}: {}\n".format(length, read_length_counts[length])

    with open("{}_bam_summary.txt".format(prefix), "w") as output:
        output.write(summary)

    return (alignments, read_length_counts)
