"""All functions that are not so useful, but still useful."""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from collections import OrderedDict
from collections import defaultdict
import csv
import errno
import itertools
import math
import os
import sys
import glob
import ntpath

from scipy import stats
import numpy as np
import pandas as pd
import six


def list_to_ranges(list_of_int):
    """Convert a list to a list of ragne object

    Parameters
    ----------

    list_of_int: list
        List of integers to be squeezed into range

    Returns
    -------

    list_of_range: list
        List of range objects


    """
    sorted_list = sorted(set(list_of_int))
    for key, group in itertools.groupby(
            enumerate(sorted_list), lambda x: x[1] - x[0]):
        group = list(group)
        yield group[0][1], group[-1][1]


def create_ideal_periodic_signal(signal_length):
    """Create ideal ribo-seq signal.

    Parameters
    ----------
    signal_length : int
                    Length of signal to create

    Returns
    -------
    signal : array_like
             1-0-0 signal

    """
    uniform_signal = np.array([4 / 6.0] * signal_length)
    uniform_signal[range(1, len(uniform_signal), 3)] = 1 / 6.0
    uniform_signal[range(2, len(uniform_signal), 3)] = 1 / 6.0
    return uniform_signal


def identify_peaks(coverage):
    """Given coverage array, find the site of maximum density"""
    return np.argmax(coverage[range(-20, -10)])


def millify(n):
    """Convert integer to human readable format.

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
    millidx = max(
        0,
        min(
            len(millnames) - 1,
            int(math.floor(0 if n == 0 else math.log10(abs(n)) / 3))))

    return '{:.1f}{}'.format(n / 10**(3 * millidx), millnames[millidx])


def mkdir_p(path):
    """Python version mkdir -p

    Parameters
    ----------

    path : str
    """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def r2(x, y):
    '''Calculate pearson correlation between two vectors.

    Parameters
    ----------
    x : array_like
        Input
    y : array_like
        Input
    '''
    return stats.pearsonr(x, y)[0]**2


def round_to_nearest(x, base=5):
    '''Round to nearest base.

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


def set_xrotation(ax, degrees):
    """Rotate labels on x-axis.

    Parameters
    ----------
    ax : matplotlib.Axes
         Axes object
    degrees : int
              Rotation degrees
    """
    for i in ax.get_xticklabels():
        i.set_rotation(degrees)


def summary_stats_two_arrays_welch(old_mean_array,
                                   new_array,
                                   old_var_array=None,
                                   old_n_counter=None,
                                   carried_forward_observations=None):
    """Average two arrays using welch's method

    Parameters
    ----------
    old_mean_array : Series
                Series of previous means with index as positions
    old_var_array : Series
                Series of previous variances with index as positions
    new_array : array like
                Series of new observations
                (Does noes
                Ciunts of number of positions at a certain index

    Returns
    -------
    m : array like
        Column wise Mean array
    var : array like
         Column wise variance

    Consider an example: [1,2,3], [1,2,3,4], [1,2,3,4,5]

    old = [1,2,3]
    new = [1,2,3,4]
    counter = [1,1,1]
    mean = [1,2,3,4] Var =[na, na, na, na], carried_fowrad = [[1,1], [2,2], [3,3], [4]]

    old = [1,2,3,4]
    new = [1,2,3,4,5]
    couter = [2,2,2,1]

    mean = [1,2,3,4,5]
    var = [0,0,0, na, na]
    carried_forward = [[], [], [], [4,4], [5]]
    """
    if not isinstance(old_mean_array, pd.Series):
        old_mean_array = pd.Series(old_mean_array)
    if not isinstance(new_array, pd.Series):
        new_array = pd.Series(new_array)
    if old_n_counter is not None and not isinstance(old_n_counter, pd.Series):
        old_n_counter = pd.Series(old_n_counter)

    len_old, len_new = len(old_mean_array), len(new_array)
    if old_n_counter is None:
        # Initlaized from current series
        old_n_counter = pd.Series(
            np.zeros(len(old_mean_array)) + 1, index=old_mean_array.index)
    if old_var_array is None:
        # Initlaized from current series
        old_var_array = pd.Series(
            np.zeros(len(old_mean_array)) + np.nan, index=old_mean_array.index)
    # Update positions counts based on new_array
    new_n_counter = old_n_counter.add(
        pd.Series(np.zeros(len(new_array)) + 1, index=new_array.index),
        fill_value=0)
    if len_old > len_new:
        len_diff = len_old - len_new
        # Pad the incoming array
        # We append NAs to the end of new_array since it will mostly be in the metagene context
        max_index = np.max(new_array.index.tolist())
        new_index = np.arange(max_index + 1, max_index + 1 + len_diff)
        new_array = new_array.append(
            pd.Series(np.zeros(len_diff) + np.nan, index=new_index),
            verify_integrity=True)
    elif len_old < len_new:
        len_diff = len_new - len_old
        # Pad the old array
        if len_old == 0:
            old_mean_array = pd.Series([])
        else:
            max_index = np.max(old_mean_array.index.tolist())
            new_index = np.arange(max_index + 1, max_index + 1 + len_diff)
            old_mean_array = old_mean_array.append(
                pd.Series(np.zeros(len_diff) + np.nan, index=new_index),
                verify_integrity=True)

    if not (old_mean_array.index == new_array.index).all():
        print('old array index: {}'.format(old_mean_array))
        print('new array index: {}'.format(new_array))
    positions_with_less_than3_obs = defaultdict(list)
    for index, counts in six.iteritems(new_n_counter):
        # Which positions has <3 counts for calculating variance
        if counts <= 3:
            # Fetch the exact observations from history
            try:
                last_observations = carried_forward_observations[index]
            except:
                # No carreid forward passed
                if not np.isnan(old_mean_array[index]):
                    last_observations = [old_mean_array[index]]
                else:
                    last_observations = []
            # Add entry from new_array only if it is not NAN
            if not np.isnan(new_array[index]):
                last_observations.append(new_array[index])
            positions_with_less_than3_obs[index] = last_observations

    # positions_with_less_than3_obs = pd.Series(positions_with_less_than3_obs)

    # delta = x_n - mean(x_{n-1})
    delta = new_array.subtract(old_mean_array)
    """
    for index, value in six.iteritems( delta ):
        if np.isnan(value):
            if not np.isnan(old_mean_array[index]):
                delta[index] = old_mean_array[index]
            else:
                delta[index] = new_array[index]
    """

    # delta = delta/n
    delta_normalized = delta.divide(new_n_counter)
    # mean(x_n) = mean(x_{n-1}) + delta/n
    new_mean_array = old_mean_array.add(delta_normalized)
    for index, value in six.iteritems(new_mean_array):
        if np.isnan(value):
            if not np.isnan(old_mean_array[index]):
                new_mean_array[index] = old_mean_array[index]
            else:
                new_mean_array[index] = new_array[index]
    #print(delta)
    #print(new_n_counter)
    #print(delta_normalized)
    #print(new_mean_array)
    # mean_difference_current = x_n - mean(x_n)
    # mean_difference_previous = x_n - mean(x_{n-1})

    mean_difference_current = new_array.fillna(0) - new_mean_array.fillna(0)
    mean_difference_previous = new_array.fillna(0) - old_mean_array.fillna(0)

    # (x_n-mean(x_n))(x_n-mean(x_{n-1})
    product = np.multiply(mean_difference_current, mean_difference_previous)

    # (n-1)S_n^2 - (n-2)S_{n-1}^2 = (x_n-mean(x_n)) (x_n-mean(x_{n-1}))

    # old_ssq = (n-1)S_{n-1}^2
    # (n-2)S_{n-1}^2
    old_sum_of_sq = (old_n_counter - 2).multiply(old_var_array.fillna(0))

    # new_ssq = (old_ssq + product)
    # (n-1) S_n^2
    new_sum_of_sq = old_sum_of_sq + product

    # if counts is less than 3, set sum of sq to NA
    new_sum_of_sq[new_n_counter < 3] = np.nan

    # if counts just became 3, compute the variance
    for index, counts in six.iteritems(new_n_counter):
        if counts == 3:
            observations = positions_with_less_than3_obs[index]
            variance = np.var(observations)
            print(index, variance)
            new_sum_of_sq[index] = variance
            # delete it from the history
            del positions_with_less_than3_obs[index]

    new_var_array = new_sum_of_sq.divide(new_n_counter - 1)
    new_var_array[new_var_array == np.inf] = np.nan
    new_var_array[new_n_counter < 3] = np.nan
    """
    for index, counts in six.iteritems(new_n_counter):
        if counts < 3:
            if not np.isnan(new_array[index]):
                if index not in list(positions_with_less_than3_obs.keys()):
                    positions_with_less_than3_obs[index] = list()
                assert index in positions_with_less_than3_obs.keys()
                positions_with_less_than3_obs[index].append(new_array[index])
    """
    return new_mean_array, new_var_array, new_n_counter, positions_with_less_than3_obs


def path_leaf(path):
    """Get path's tail from a filepath"""
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)


def parse_star_logs(f):
    """Parse star logs into a dict

    Parameters
    ----------
    f : str
        Path to starlogs.final.out file

    Returns
    -------
    star_info : dict
                Dict with necessary records parsed
    """
    ANNOTATIONS = [
        'Total reads', 'Uniquely Mapped', 'Uniquely Mapped %',
        'Multi Mapped %', 'Unmapped %', 'Multi Mapped'
    ]
    star_info = OrderedDict()
    with open(f) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('Number of input reads'):
                star_info[ANNOTATIONS[0]] = int(line.strip().split('\t')[1])
            elif line.startswith('Uniquely mapped reads number'):
                star_info[ANNOTATIONS[1]] = int(line.strip().split('\t')[1])
            elif line.startswith('Uniquely mapped reads %'):
                star_info[ANNOTATIONS[2]] = round(
                    float(line.strip('%').split('\t')[1]), 2)
            elif line.startswith('Number of reads mapped to multiple loci'):
                star_info[ANNOTATIONS[5]] = int(line.strip().split('\t')[1])
            elif line.startswith('Number of reads mapped to too many loci'):
                star_info[ANNOTATIONS[5]] += int(line.strip().split('\t')[1])
            elif line.startswith('% of reads mapped to multiple loci'):
                star_info[ANNOTATIONS[3]] = round(
                    float(line.strip('%').split('\t')[1]), 2)
            elif line.startswith('% of reads mapped to too many loci'):
                star_info[ANNOTATIONS[3]] += round(
                    float(line.strip('%').split('\t')[1]), 2)
            elif line.startswith('% of reads unmapped: too many mismatches'):
                star_info[ANNOTATIONS[4]] = round(
                    float(line.strip('%').split('\t')[1]), 2)
            elif line.startswith('% of reads unmapped: too short'):
                star_info[ANNOTATIONS[4]] += round(
                    float(line.strip('%').split('\t')[1]), 2)
            elif line.startswith('% of reads unmapped: other'):
                star_info[ANNOTATIONS[4]] += round(
                    float(line.strip('%').split('\t')[1]), 2)

    star_info = {
        key: round(star_info[key], 2)
        for key in list(star_info.keys())
    }
    return star_info


def get_strandedness(filepath):
    """Parse output of infer_experiment.py from RSeqC to get strandedness.

    Parameters
    ----------
    filepath : str
               Path to infer_experiment.py output

    Returns
    -------
    strandedness : str
                   reverse or forward or none
    """
    with open(filepath) as f:
        data = f.read()
    splitted = [x.strip() for x in data.split('\n') if len(x.strip()) >= 1]
    strandedness = None
    assert splitted[0] == 'This is SingleEnd Data'
    few_percentage = None
    rev_percentage = None
    for line in splitted[1:]:
        if 'Fraction of reads failed to determine:' in line:
            continue
        elif 'Fraction of reads explained by "++,--":' in line:
            fwd_percentage = float(line.split(':')[1])
        elif 'Fraction of reads explained by "+-,-+":' in line:
            rev_percentage = float(line.split(':')[1])

    assert rev_percentage is not None
    assert fwd_percentage is not None

    ratio = fwd_percentage / rev_percentage

    if np.isclose([ratio], [1]):
        return 'none'
    elif ratio >= 0.5:
        return 'forward'
    else:
        return 'reverse'
