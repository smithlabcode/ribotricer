from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
import pandas as pd
from scipy.stats.mstats import ks_2samp
import statsmodels.api as sm


def KDE(values):
    """Perform Univariate Kernel Density Estimation.

    Wrapper utility around statsmodels for quick KDE
    TODO: scikit-learn has a faster implementation (?)

    Parameters
    ----------
    values : array like


    Returns
    -------

    support : array_like
    cdf : array_like

    """
    density = sm.nonparametric.KDEUnivariate(values)
    density.fit()
    return density.support, density.cdf


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


def series_cdf(series):
    """Calculate cdf of series preserving the index

    Parameters
    ----------
    series : series like

    Returns
    -------
    cdf : series

    """
    index = series.index.tolist()
    total = series.sum(skipna=True)
    cdf = series.cumsum(skipna=True) / total
    return cdf


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
    if not isinstance(a, pd.Series):
        a = pd.Series(a)
    if not isinstance(b, pd.Series):
        b = pd.Series(b)

    cdf_a = series_cdf(a)
    cdf_b = series_cdf(b)
    effect_size, p = ks_2samp(a, b, alternative='greater')
    D = cdf_a.subtract(cdf_b).abs().idxmax()
    return D, effect_size, p, cdf_a, cdf_b
