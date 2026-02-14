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

from __future__ import annotations

import ntpath
import pathlib
import sys
from typing import TYPE_CHECKING

from .interval import Interval

if TYPE_CHECKING:
    import pysam

# Source: https://broadinstitute.github.io/picard/explain-flags.html
__SAM_NOT_UNIQ_FLAGS__: list[int] = [4, 20, 256, 272, 2048]


def is_read_uniq_mapping(read: pysam.AlignedSegment) -> bool | None:
    """Check if read is uniquely mappable.

    Parameters
    ----------
    read : pysam.AlignedSegment
        A pysam alignment object.

    Returns
    -------
    bool | None
        True if uniquely mapping, False if not, None if unable to determine.

    Notes
    -----
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
        elif read.mapping_quality < 1 or read.flag in __SAM_NOT_UNIQ_FLAGS__:
            return False
        else:
            sys.stdout.write(
                "WARNING: ribotricer was unable to detect any tags for "
                "determining multimapping status. All the reads will be "
                "treated as uniquely mapping\n"
            )
            return None


def merge_intervals(intervals: list[Interval]) -> list[Interval]:
    """Merge overlapping intervals.

    Parameters
    ----------
    intervals : list[Interval]
        List of intervals to merge.

    Returns
    -------
    list[Interval]
        Sorted and merged intervals.
    """
    intervals = sorted(intervals, key=lambda x: x.start)
    merged_intervals: list[Interval] = []
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


def mkdir_p(path: str) -> None:
    """Make directory even if it exists.

    Parameters
    ----------
    path : str
        Path to directory to create.
    """
    pathlib.Path(path).mkdir(parents=True, exist_ok=True)


def path_leaf(path: str) -> str:
    """Get path's tail from a filepath.

    Parameters
    ----------
    path : str
        File path.

    Returns
    -------
    str
        The tail (filename) portion of the path.
    """
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)


def parent_dir(path: str) -> str:
    """Get path's parent directory from a filepath.

    Parameters
    ----------
    path : str
        File path.

    Returns
    -------
    str
        The parent directory portion of the path.
    """
    head, tail = ntpath.split(path)
    return head


def _clean_input(comma_string: str) -> list[str]:
    """Clean comma separated option inputs in CLI.

    Parameters
    ----------
    comma_string : str
        Comma-separated string of values.

    Returns
    -------
    list[str]
        List of stripped string values.
    """
    return [term.strip(" ") for term in comma_string.split(",")]


def collapse_coverage_to_codon(coverage: list[int]) -> list[int]:
    """Collapse nucleotide level coverage to codon level.

    Parameters
    ----------
    coverage : list[int]
        Nucleotide level counts.

    Returns
    -------
    list[int]
        Coverage collapsed to codon level.
    """
    codon_coverage = [
        sum(coverage[current : current + 3]) for current in range(0, len(coverage), 3)
    ]
    return codon_coverage
