"""Statistics related functions
"""
import numpy as np
import scipy as sp
from .interval import Interval
from scipy import stats
from scipy import signal
import scipy.integrate as integrate
from math import sin, cos, pi, sqrt


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
    """Calculate coherence with an idea ribo-seq signal

    Parameters
    ----------
    values : array like
             List of value

    Returns
    -------
    coh : float
          Periodicity score calculated as
          coherence between input and idea 1-0-0 signal

    pval: float
          p value for the coherence

    valid: int
           number of valid codons used for calculation

    """
    coh, pval, valid = 0.0, 1.0, -1
    for frame in [0, 1, 2]:
        values = original_values[frame:]
        normalized_values = []
        i = 0
        while i + 2 < len(values):
            if values[i] == values[i + 1] == values[i + 2] == 0:
                i += 3
                continue
            real = values[i] + values[i + 1] * cos(
                2 * pi / 3) + values[i + 2] * cos(4 * pi / 3)
            image = values[i + 1] * sin(2 * pi / 3) + values[i + 2] * sin(
                4 * pi / 3)
            norm = sqrt(real**2 + image**2)
            if norm == 0:
                norm = 1
            normalized_values += [
                values[i] / norm, values[i + 1] / norm, values[i + 2] / norm
            ]
            i += 3

        length = len(normalized_values) // 3 * 3
        if length == 0:
            return (0.0, 1.0, 0)
        normalized_values = normalized_values[:length]
        uniform_signal = [1, 0, 0] * (len(normalized_values) // 3)
        f, Cxy = signal.coherence(
            normalized_values,
            uniform_signal,
            window=[1.0, 1.0, 1.0],
            nperseg=3,
            noverlap=0)
        try:
            periodicity_score = Cxy[np.argwhere(np.isclose(f, 1 / 3.0))[0]][0]
            periodicity_pval = pvalue(periodicity_score, length // 3)
        except:
            periodicity_score = 0.0
            periodicity_pval = 1.0
        if periodicity_score > coh:
            coh = periodicity_score
            pval = periodicity_pval
            valid = length
        if valid == -1:
            valid = length
    return coh, pval, valid
