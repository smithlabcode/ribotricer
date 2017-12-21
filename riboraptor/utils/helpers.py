"""All functions that are not so useful, but still useful."""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import math
from scipy import stats
import numpy as np


def millify(n):
    """Convert integer to human readable format

    Parameters
    ----------
    n : int

    Returns
    -------
    millidx : str
              Formatted integer
    """
    millnames = ['', ' K', ' M', ' B', ' T']
    # Source: http://stackoverflow.com/a/3155023/756986
    n = float(n)
    millidx = max(0, min(len(millnames) - 1,
                         int(math.floor(0 if n == 0 else math.log10(abs(n)) / 3))))

    return '{:.0f}{}'.format(n / 10**(3 * millidx), millnames[millidx])


def create_ideal_periodic_signal(signal_length):
    """Create ideal ribo-seq signal

    Parameters
    ---------
    signal_length : int

    Returns
    -------

    signal : array_like
             1-0-0 signal

    """
    uniform_signal = np.array([4 / 6.0] * signal_length)
    uniform_signal[range(1, len(uniform_signal), 3)] = 1 / 6.0
    uniform_signal[range(2, len(uniform_signal), 3)] = 1 / 6.0
    return uniform_signal


def set_rotation(ax, degrees):
    """Rotate labels on axis.
    Parameters
    ----------
    ax : matplotlib.Axes
    degrees : int
              Rotation degrees
    """
    labels = ax.get_xticklabels()
    for i in labels:
        i.set_rotation(degrees)


def round_to_nearest(x, base=5):
    '''Round to nearest base

    Parameters
    ----------
    x : float
        Input

    Returns
    -------
    v : int
        Output
    '''
    return int(base * round(float(x) / base))


def r2(x, y):
    '''Calculate pearson correlation between two vectors

    Parameters
    ----------
    x : array_like
        Input
    y : array_like
        Input
    '''
    return stats.pearsonr(x, y)[0] ** 2


def identify_peaks(coverage):
    """Given coverage array, find the site of maximum density"""
    return np.argmax(coverage[range(-30, 0)])