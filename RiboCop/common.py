"""Utilities for common usage
"""
import numpy as np
import scipy as sp
from .interval import Interval
from .test_func import wilcoxon_greater, combine_pvals
from scipy import stats
from scipy import signal
from math import sin, cos, pi, sqrt

def phase_vector_pdf(x, k, form='scipy'):
    """Get PDF for distribution of uniform random walk

    Parameters
    ----------
    x: array_like
       An array with values of $x$ at which P(X) is to
       be calculated. Example: np.linspace(1, 100, 10000)
    k: int
       Number of vectors projected (total number of segments/windows)
       It is equal to $N$ in our coherence derivation

    Returns
    -------
    p: array_like
       Array with PDF calcualted at point of $x$

    Reference: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1096409&tag=1
    """
    if k == 1:
        return [0]*len(x)
    if k == 2:
        p = (1/np.pi) * (1/np.sqrt(np.multiply(x, 4-x)))
        # Set p to zero if it is outside the range of [0, 4]
        p = np.multiply(( p <= 4) & (p >= 0), p)
    elif k >= 3:
        # Warning this is currently not stable and will require
        # recursion
        if form != 'scipy':
            p = 1/(k-1) * np.multiply(np.exp(-(x+1)/(k-1)) , sp.special.jv(0, 2*np.sqrt(x)/(k-1)))
            p = np.multiply((p >= 0), p)
        else:
            x = 2 * x /(k - 1)
            df = 2
            nc = 2 / (k -1)
            p = 2 * stats.ncx2.pdf(x, df, nc)/(k-1)
    return p

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


def pvalue(x, N):
    """Calculate p-value for coherence score

    Parameters
    ----------
    x: double
       coherence score
    N: int
       number of valid codons

    Returns
    -------
    pval: double
          p-value for the coherence score
    """
    df, nc = 2, 2.0 / (N - 1)
    x = 2 * N**2 * x / (N - 1)
    return stats.ncx2.sf(x, df, nc)


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
    coh, pval, valid = 0.0, 1.0, -1
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
            i += 3

        length = len(normalized_values) // 3 * 3
        if length == 0:
            return (0.0, 1.0, 0)
        normalized_values = normalized_values[:length]
        uniform_signal = [1, 0, 0] * (len(normalized_values) // 3)
        f, Cxy = signal.coherence(
            normalized_values, uniform_signal, window=[1.0,1.0,1.0], nperseg=3, noverlap=0)
        try:
            periodicity_score = Cxy[np.argwhere(np.isclose(f, 1 / 3.0))[0]][0]
            periodicity_pval = pvalue(periodicity_score, length // 3)
        except:
            periodicity_score = 0.0
            periodicity_pval = 1.0
        if periodicity_pval < pval:
            coh = periodicity_score
            pval = periodicity_pval
            valid = length
        if valid == -1:
            valid = length
    return coh, pval, valid


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
