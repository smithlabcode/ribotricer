"""Utilities for common usage
"""
import numpy as np
from .interval import Interval
from .test_func import wilcoxon_greater, combine_pvals
from scipy import stats
from scipy import signal


def cal_periodicity(values):
    coh = coherence(values)
    pval = wilcoxon(values)
    return (coh, pval)


def wilcoxon(values):
    length = len(values) // 3 * 3
    values = values[:length]
    f0 = values[0:length:3]
    f1 = values[1:length:3]
    f2 = values[2:length:3]
    final_pv1=final_pv2=final_pv=1.0

    pv1 = wilcoxon_greater(f0, f1)
    pv2 = wilcoxon_greater(f0, f2)
    pv = combine_pvals(np.array([pv1.pvalue, pv2.pvalue]))
    if pv < final_pv:
        final_pv = pv
        final_pv1 = pv1
        final_pv2 = pv2

    pv1 = wilcoxon_greater(f1, f0)
    pv2 = wilcoxon_greater(f1, f2)
    pv = combine_pvals(np.array([pv1.pvalue, pv2.pvalue]))
    if pv < final_pv:
        final_pv = pv
        final_pv1 = pv1
        final_pv2 = pv2

    pv1 = wilcoxon_greater(f2, f0)
    pv2 = wilcoxon_greater(f2, f1)
    pv = combine_pvals(np.array([pv1.pvalue, pv2.pvalue]))
    if pv < final_pv:
        final_pv = pv
        final_pv1 = pv1
        final_pv2 = pv2
    return final_pv

def coherence(values):
    """Calculate coherence and an idea ribo-seq signal

    Parameters
    ----------
    values : array like
             List of values

    Returns
    -------
    periodicity : float
                  Periodicity score calculated as
                  coherence between input and idea 1-0-0 signal

    f: array like
       List of frequencies

    Cxy: array like
         List of coherence at the above frequencies

    """
    length = len(values) // 3 * 3
    values = values[:length]
    uniform_signal = [0.7, 0.2, 0.1] * (length // 3)
    mean_centered_values = values - np.nanmean(values)
    normalized_values = mean_centered_values / \
        np.max(np.abs(mean_centered_values))

    mean_centered_values = uniform_signal - np.nanmean(uniform_signal)
    uniform_signal = mean_centered_values / \
        np.max(np.abs(uniform_signal))
    f, Cxy = signal.coherence(
        normalized_values, uniform_signal, nperseg=30, noverlap=27)
    periodicity_score = Cxy[np.argwhere(np.isclose(f, 1 / 3.0))[0]][0]
    # return periodicity_score, f, Cxy
    return periodicity


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
