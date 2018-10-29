"""Utilities for common usage
"""
import numpy as np
from .interval import Interval
from .test_func import wilcoxon_greater, combine_pvals
from scipy import stats
from scipy import signal
from math import sin, cos, pi, sqrt


def cal_periodicity(values):
    coh, nonzero = coherence(values)
    # pval, nonzero = wilcoxon(values)
    return (coh, 1.0, nonzero)


def wilcoxon(values):
    length = len(values) // 3 * 3
    values = values[:length]
    f0 = np.array(values[0:length:3])
    f1 = np.array(values[1:length:3])
    f2 = np.array(values[2:length:3])
    final_pv1 = final_pv2 = final_pv = 1.0
    nonzero = 0

    pv1 = wilcoxon_greater(f0, f1)
    pv2 = wilcoxon_greater(f0, f2)
    pv = combine_pvals(np.array([pv1.pvalue, pv2.pvalue]))
    if pv < final_pv:
        final_pv = pv
        final_pv1 = pv1
        final_pv2 = pv2
        nonzero = np.flatnonzero(f0).size

    pv1 = wilcoxon_greater(f1, f0)
    pv2 = wilcoxon_greater(f1, f2)
    pv = combine_pvals(np.array([pv1.pvalue, pv2.pvalue]))
    if pv < final_pv:
        final_pv = pv
        final_pv1 = pv1
        final_pv2 = pv2
        nonzero = np.flatnonzero(f1).size

    pv1 = wilcoxon_greater(f2, f0)
    pv2 = wilcoxon_greater(f2, f1)
    pv = combine_pvals(np.array([pv1.pvalue, pv2.pvalue]))
    if pv < final_pv:
        final_pv = pv
        final_pv1 = pv1
        final_pv2 = pv2
        nonzero = np.flatnonzero(f2).size
    return final_pv, nonzero


def repeat_codon(x, times):
    ans = []
    for i in range(0, len(x), 3):
        ans += (x[i:i+3] * times)
    return ans


def coherence(original_values):
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
    p, valid, predict_frame = 0.0, -1, 0
    for frame in [0, 1, 2]:
        values = original_values[frame:]
        normalized_values = []
        i = 0
        while i + 2 < len(values):
            if values[i] == values[i+1] == values[i+2] == 0:
                i += 3
                continue
            real = values[i] + values[i+1] * cos(2*pi/3) + values[i+2] * cos(4*pi/3)
            image = values[i+1] * sin(2*pi/3) + values[i+2] * sin(4*pi/3)
            norm = sqrt(real ** 2 + image ** 2)
            if norm == 0:
                norm = 1
            normalized_values += [values[i] / norm, values[i+1] / norm, values[i+2] / norm]
            # if values[i] == 0:
            #     normalized_values += [values[i], values[i+1], values[i+2]]
            # else:
            #     normalized_values += [1.0, values[i+1] / values[i], values[i+2] / values[i]]
            i += 3

        length = len(normalized_values) // 3 * 3
        if length == 0:
            return (-1, -1, -1)
        normalized_values = normalized_values[:length]
        uniform_signal = [1, 0, 0] * (len(normalized_values) // 3)
        f, Cxy = signal.coherence(
            normalized_values, uniform_signal, window=[1.0,1.0,1.0], nperseg=3, noverlap=0)
        try:
            periodicity_score = Cxy[np.argwhere(np.isclose(f, 1 / 3.0))[0]][0]
        except:
            periodicity_score = 0.0
        if periodicity_score > p:
            p = periodicity_score
            valid = length
            predict_frame = frame
        if valid == -1:
            valid = length
    return p, valid, predict_frame


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
