from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
import pandas as pd
from scipy.stats.mstats import ks_2samp


def calculate_cdf(data):
    """Calculate CDF given data points

    Parameters
    ----------
    data : array-like
        Input values

    Returns
    -------
    cdf : series
        Cumulative distribution funvtion calculated at indexed points

    """
    data = pd.Series(data)
    data = data.fillna(0)
    total = np.nansum(data)
    index = data.index.tolist()
    cdf = []
    for i in range(min(index), max(index)):
        cdf.append(np.sum(data[range(min(index), i)]) / total)
    return pd.Series(cdf, index=index[1:])


def KS_test(a, b):
    """Perform KS test between a and b values

    Parameters
    ----------
    a, b : array-like
           Input

    Returns
    -------
    D : int
        KS D statistic
    effect_size : float
                  maximum difference at point of D-statistic
    cdf_a, cdf_b : float
                   CDF of a, b

    Note:    By default this method does testing for alternative=lesser implying
    that the test will reject H0 when the CDf of b is 'above' a


    """
    cdf_a = calculate_cdf(a)
    cdf_b = calculate_cdf(b)
    effect_size, p = ks_2samp(a, b, alternative='greater')
    D = np.argmax(np.abs(cdf_a - cdf_b)) - 1
    # a.index[np.argmax(np.array(cdf_a)-np.array(cdf_b))]
    return D, effect_size, p, cdf_a, cdf_b
