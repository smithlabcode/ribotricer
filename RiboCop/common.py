"""Utilities for common usage
"""
import numpy as np
from .interval import Interval
from scipy import stats


def cal_periodicity(values):
    repeats = len(values) // 3
    total = repeats * 3
    values = values[:total]
    corr, pval = 0, 1.0
    frame0 = np.array([1, 0, 0] * repeats)
    r, p = stats.pearsonr(values, frame0)
    if abs(r) > corr:
        corr, pval = abs(r), p

    frame1 = np.array([0, 1, 0] * repeats)
    r, p = stats.pearsonr(values, frame1)
    if abs(r) > corr:
        corr, pval = abs(r), p

    frame2 = np.array([0, 0, 1] * repeats)
    r, p = stats.pearsonr(values, frame2)
    if abs(r) > corr:
        corr, pval = abs(r), p
    return (corr, pval)


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
        nh_count = tags['NH']
    except KeyError:
        # Reliable in case of STAR
        if read.mapping_quality == 255:
            return True
        if read.mapping_quality < 1:
            return False
        # NH tag not set so rely on flags
        if read.flag in __SAM_NOT_UNIQ_FLAGS__:
            return False
        else:
            raise RuntimeError('Malformed BAM?')
    if nh_count == 1:
        return True
    return False


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
        to_merge = Interval(intervals[i].chrom, intervals[i].start,
                            intervals[i].end, intervals[i].strand)
        while (i + 1 < len(intervals)
               and intervals[i + 1].start <= to_merge.end):
            to_merge.end = max(to_merge.end, intervals[i + 1].end)
            i += 1
        merged_intervals.append(to_merge)
        i += 1
    return merged_intervals
