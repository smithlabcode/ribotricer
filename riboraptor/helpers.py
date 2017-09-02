from scipy import stats

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

