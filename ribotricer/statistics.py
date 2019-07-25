"""Statistics related functions
"""
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

from math import sin, cos, pi, sqrt
import numpy as np
from scipy import stats
from scipy import signal


def phasescore(profile):
    """Calculate phase score for a RPF profile

    Parameters
    ----------
    profile: array like
             list of value

    Returns
    -------
    phase_score: double
                 phase score for RPF profile
    nonempty: int
              numbere of non-empty codons
    """
    phase_score, nonempty = 0.0, -1
    for frame in [0, 1, 2]:
        values = profile[frame:]
        cur_phase_score = 0.0
        x = y = 0.0
        valid = 0
        i = 0
        while i + 2 < len(values):
            if values[i] == values[i + 1] == values[i + 2] == 0:
                i += 3
            else:
                valid += 1
                cur_x = (
                    values[i] * 1
                    + values[i + 1] * cos(-2 * pi / 3)
                    + values[i + 2] * cos(-4 * pi / 3)
                )
                cur_y = (
                    values[i] * 0
                    + values[i + 1] * sin(-2 * pi / 3)
                    + values[i + 2] * sin(-4 * pi / 3)
                )
                norm = sqrt(cur_x * cur_x + cur_y * cur_y)
                if norm > 0:
                    x += cur_x / norm
                    y += cur_y / norm
                i += 3
        if valid > 0:
            cur_phase_score = sqrt(x * x + y * y) / valid
        if nonempty == -1:
            nonempty = valid
        if cur_phase_score > phase_score:
            phase_score = cur_phase_score
            nonempty = valid
    return phase_score, nonempty


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
    x = 2 * N ** 2 * x / (N - 1)
    return stats.ncx2.sf(x, df, nc)


def coherence(original_values):
    """Calculate coherence with an idea ribo-seq signal,
    this function is equivalent to phasescore implemented
    above, but we expect scipy has better handling of numerical
    precision.

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
