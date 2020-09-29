"""Statistics related functions
"""
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

from math import sin, cos, pi, sqrt
import warnings

import numpy as np
from scipy import stats
from scipy import signal


def pvalue(x, N):
    """Calculate p-value for phase score

    Parameters
    ----------
    x: double
       phase score
    N: int
       number of valid codons

    Returns
    -------
    pval: double
          p-value for the phase score
    """
    df, nc = 2, 2.0 / (N - 1)
    x = 2 * N ** 2 * x / (N - 1)
    return stats.ncx2.sf(x, df, nc)


def phasescore(original_values):
    """Calculate phase score of a given signal.

    Parameters
    ----------
    values : array like
             List of value

    Returns
    -------
    coh : float
          Periodicity score calculated as
          coherence between input and idea 1-0-0 signal

    valid: int
           number of valid codons used for calculation

    """
    coh, valid = 0.0, -1
    for frame in [0, 1, 2]:
        values = original_values[frame:]
        normalized_values = []
        i = 0
        while i + 2 < len(values):
            if values[i] == values[i + 1] == values[i + 2] == 0:
                i += 3
            else:
                real = (
                    values[i]
                    + values[i + 1] * cos(2 * pi / 3)
                    + values[i + 2] * cos(4 * pi / 3)
                )
                image = values[i + 1] * sin(2 * pi / 3) + values[i + 2] * sin(
                    4 * pi / 3
                )
                norm = sqrt(real ** 2 + image ** 2)
                if norm == 0:
                    norm = 1
                normalized_values += [
                    values[i] / norm,
                    values[i + 1] / norm,
                    values[i + 2] / norm,
                ]
                i += 3

        length = len(normalized_values) // 3 * 3
        if length == 0:
            coh, valid = (0.0, 0)
        else:
            normalized_values = normalized_values[:length]
            uniform_signal = [1, 0, 0] * (len(normalized_values) // 3)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                f, Cxy = signal.coherence(
                    normalized_values,
                    uniform_signal,
                    window=[1.0, 1.0, 1.0],
                    nperseg=3,
                    noverlap=0,
                )
                periodicity_score = Cxy[np.argwhere(np.isclose(f, 1 / 3.0))[0]][0]
                if periodicity_score > coh:
                    coh = periodicity_score
                    valid = length // 3
                if valid == -1:
                    valid = length // 3
    return np.sqrt(coh), valid
