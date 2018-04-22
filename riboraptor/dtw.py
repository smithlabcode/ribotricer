from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
import matplotlib.pyplot as plt

from scipy.spatial.distance import cdist
from matplotlib.ticker import NullFormatter


def dtw(X, Y, metric='euclidean', ddtw=False, ddtw_order=1):
    """

    Parameters
    ----------
    X : array_like
        M x D matrix
    Y : array_like
        N x D matrix
    metric : string
             The distance metric to use.
             Can be :
             'braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation',
             'cosine', 'dice', 'euclidean', 'hamming', 'jaccard', 'kulsinski',
             'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao',
             'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean',
             'wminkowski', 'yule'.
             See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.cdist.html
    ddtw : bool
           Should use derivative DTW where the distance matrix is created
           using the derivate values at each point rather than the point themselves
    ddtw_order : int [1,2]
                 First order uses one difference method
                 Second order uses np.gradient for an approximation upto second order
    Returns
    -------
    total_cost : float
                 Total (minimum) cost of warping
    pointwise_cost :  array_like
                      M x N matrix with cost at each (i, j)
    accumulated_cost : array_like
                       M x N matrix with (minimum) cost accumulated till (i,j)
                       having started from (0, 0)

    """
    X = list(X)
    Y = list(Y)
    X = np.array(X)
    Y = np.array(Y)
    if len(X.shape) == 1:
        # Reshape to N x 1 form
        X = X[:, np.newaxis]
    if len(Y.shape) == 1:
        Y = Y[:, np.newaxis]
    # m = X length
    # n = Y length
    m, n = X.shape[0], Y.shape[0]
    assert X.shape[1] == Y.shape[1]
    d = X.shape[1]
    D = np.zeros((m + 1, n + 1))
    D[1:, 0] = np.inf
    D[0, 1:] = np.inf
    if not ddtw:
        D[1:, 1:] = cdist(X, Y, metric)
    else:
        if ddtw_order == 1:
            X_mod = np.zeros(X.shape)
            Y_mod = np.zeros(Y.shape)
            # shape will be M-1 x D
            X_diff = np.diff(X, axis=0)
            # shape will be M-1 x D
            Y_diff = np.diff(Y, axis=0)

            X_mod[1:m - 1, :] = 0.5 * (
                X_diff[:m - 2, :] + 0.5 * X_diff[1:m - 1])

            Y_mod[1:n - 1, :] = 0.5 * (
                Y_diff[:n - 2, :] + 0.5 * Y_diff[1:n - 1])

            X_mod[0, :] = X_mod[1, :]
            X_mod[m - 1, :] = X_mod[m - 2, :]

            Y_mod[0, :] = Y_mod[1, :]
            Y_mod[n - 1, :] = Y_mod[n - 2, :]
        elif ddtw_order == 2:
            X_mod = np.gradient(X, axis=0)
            Y_mod = np.gradient(Y, axis=0)
        else:
            raise NotImplemented('Not implemented order {}'.format(ddtw_order))
        D[1:, 1:] = cdist(X_mod, Y_mod, metric)
    pointwise_cost = D[1:, 1:].copy()
    for i in range(0, m):
        for j in range(0, n):
            cost = D[i + 1, j + 1]
            D[i + 1, j + 1] = cost + min(D[i, j + 1], D[i + 1, j], D[i, j])
    accumulated_cost = D[1:, 1:]
    total_cost = D[m, n] / sum(D.shape)
    return total_cost, pointwise_cost, accumulated_cost


def get_path(D):
    """Traceback path of minimum cost

    Given accumulated cost matrix D,
    trace back the minimum cost path

    Parameters
    -----------

    D : array_like
        M x N matrix as obtained from `accumulated_cost` using:
        total_cost, pointwise_cost, accumulated_cost = dtw(X, Y, metric='euclidean')

    Returns
    -------
    traceback_x, traceback_x : array_like
                               M x 1 and N x 1 array containing  indices of movement
                               starting from (0, 0) going to (M-1, N-1)
    """
    m, n = D.shape
    print(m, n)
    m = m - 1
    n = n - 1
    # Starting point is the end point
    traceback_x, traceback_y = [m], [n]
    while (m > 0 and n > 0):
        min_idx = np.argmin([D[m - 1, n - 1], D[m, n - 1], D[m - 1, n]])
        if min_idx == 0:
            # move diagonally
            m = m - 1
            n = n - 1
        elif min_idx == 1:
            # move vertically
            n = n - 1
        else:
            # move horizontally
            m = m - 1
        traceback_x.insert(0, m)
        traceback_y.insert(0, n)
    # End point is the starting point
    traceback_x.insert(0, 0)
    traceback_y.insert(0, 0)
    print(len(traceback_x), len(traceback_y))
    return np.array(traceback_x), np.array(traceback_y)


def plot_warped_timeseries(x,
                           y,
                           pointwise_cost,
                           accumulated_cost,
                           path,
                           colormap=plt.cm.Blues,
                           linecolor=CBB_PALETTE[-2]):
    nullfmt = NullFormatter()
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02
    rect_heatmap = [left, bottom, width, height]
    rect_x = [left, bottom_h, width, 0.2]
    rect_y = [left_h, bottom, 0.2, height]

    fig = plt.figure(1, figsize=(8, 8))

    axHeatmap = plt.axes(rect_heatmap)
    axX = plt.axes(rect_x, sharex=axHeatmap)
    axY = plt.axes(rect_y, sharey=axHeatmap)

    # no labels
    axX.xaxis.set_major_formatter(nullfmt)
    axY.yaxis.set_major_formatter(nullfmt)

    axY.plot(list(y), range(0, len(y)), color=CBB_PALETTE[2])
    axX.plot(list(x), color=CBB_PALETTE[1])
    axHeatmap.imshow(
        accumulated_cost.T,
        origin='lower',
        cmap=colormap,
        interpolation='nearest')
    axHeatmap.plot(path[0], path[1], '-x', color=linecolor)
    #axHeatmap.xlim((-0.5, accumulated_cost.shape[0]-0.5))
    #axHeatmap.ylim((-0.5, accumulated_cost.shape[1]-0.5))
    return fig
