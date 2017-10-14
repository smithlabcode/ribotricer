from scipy import stats
import numpy as np


def create_uniform_signal(signal_length):
    uniform_signal = np.array([4 / 6.0] * signal_length)
    uniform_signal[range(1, len(uniform_signal), 3)] = 1 / 6.0
    uniform_signal[range(2, len(uniform_signal), 3)] = 1 / 6.0
    return uniform_signal


def set_rotation(ax):
    labels = ax.get_xticklabels()
    for i in labels:
        i.set_rotation(45)


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
    return int(base * round(float(x)/base))


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
