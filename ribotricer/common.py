"""Utilities for common usage"""
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

import ntpath
import pathlib
import sys
from .interval import Interval

# Source: https://broadinstitute.github.io/picard/explain-flags.html
__SAM_NOT_UNIQ_FLAGS__ = [4, 20, 256, 272, 2048]


def is_read_uniq_mapping(read):
    """Check if read is uniquely mappable.

    Parameters
    ----------
    read : pysam.Alignment.fetch object


    Most reliable: ['NH'] tag
    """
    # Filter out secondary alignments
    if read.is_secondary:
        return False
    tags = dict(read.get_tags())
    try:
        nh_count = tags["NH"]
        return nh_count == 1
    except KeyError:
        # Reliable in case of STAR
        if read.mapping_quality == 255:
            return True
        elif read.mapping_quality < 1:
            return False
        # NH tag not set so rely on flags
        elif read.flag in __SAM_NOT_UNIQ_FLAGS__:
            return False
        else:
            sys.stdout.write(
                "WARNING: ribotricer was unable to detect any tags for determining multimapping status. All the reads will be treated as uniquely mapping\n"
            )


def merge_intervals(intervals):
    """
    Parameters
    ----------
    intervals: List[Interval]

    Returns
    -------
    merged_intervals: List[Interval]
                      sorted and merged intervals
    """

    intervals = sorted(intervals, key=lambda x: x.start)
    merged_intervals = []
    i = 0
    while i < len(intervals):
        to_merge = Interval(
            intervals[i].chrom,
            intervals[i].start,
            intervals[i].end,
            intervals[i].strand,
        )
        while i + 1 < len(intervals) and intervals[i + 1].start <= to_merge.end:
            to_merge.end = max(to_merge.end, intervals[i + 1].end)
            i += 1
        merged_intervals.append(to_merge)
        i += 1
    return merged_intervals


def mkdir_p(path):
    """Make directory even if it exists.

    Parameters
    ----------
    path: str
    """
    pathlib.Path(path).mkdir(parents=True, exist_ok=True)


def path_leaf(path):
    """Get path's tail from a filepath"""
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)


def parent_dir(path):
    """Get path's tail from a filepath"""
    head, tail = ntpath.split(path)
    return head


def _clean_input(comma_string):
    """Clean comma separated option inputs in CLI"""
    return list(map(lambda term: term.strip(" "), comma_string.split(",")))


def collapse_coverage_to_codon(coverage):
    """Collapse nucleotide level coverage to codon level.

    Parameters
    ----------
    coverage: list
              Nucleotide level counts
    Returns
    -------
    codon_coverage: list
                    Coverage collapsed to codon level
    """
    codon_coverage = [
        sum(coverage[current : current + 3]) for current in range(0, len(coverage), 3)
    ]
    return codon_coverage
