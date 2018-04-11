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
import pickle

from scipy import stats
import numpy as np
import pandas as pd
import six

CBB_PALETTE = [
    "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7"
]


def _fix_bed_coltype(bed):
    """Fix bed chrom and name columns to be string

    This is necessary since the chromosome numbers are often interpreted as int
    """
    bed['chrom'] = bed['chrom'].astype(str)
    bed['name'] = bed['name'].astype(str)
    return bed


def check_file_exists(filepath):
    """Check if file exists.

    Parameters
    ----------
    filepath : str
               Path to file
    """
    if os.path.isfile(os.path.abspath(filepath)):
        return True
    return False


def list_to_ranges(list_of_int):
    """Convert a list to a list of range object

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
    uniform_signal[list(range(1, len(uniform_signal), 3))] = 1 / 6.0
    uniform_signal[list(range(2, len(uniform_signal), 3))] = 1 / 6.0
    return uniform_signal


def identify_peaks(coverage):
    """Given coverage array, find the site of maximum density"""
    return np.argmax(coverage[list(range(-18, -10))])


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


def parse_star_logs(infile, outfile=None):
    """Parse star logs into a dict

    Parameters
    ----------
    infile : str
        Path to starlogs.final.out file

    Returns
    -------
    star_info : dict
                Dict with necessary records parsed
    """
    ANNOTATIONS = [
        'total_reads', 'uniquely_mapped', 'uniquely_mapped_percent',
        'multi_mapped_percent', 'unmapped_percent', 'multi_mapped'
    ]
    star_info = OrderedDict()
    with open(infile) as fh:
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
    if outfile is None:
        return star_info
    filename = path_leaf(infile)
    filename = filename.strip('Log.final.out')
    counts_df = pd.DataFrame.from_dict(star_info, orient='index').T
    counts_df.index = [filename]
    counts_df.to_csv(outfile, sep=str('\t'), index=True, header=True)
    return counts_df


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


def load_pickle(filepath):
    """Read pickled files easy in Python 2/3"""
    if sys.version_info > (3, 0):
        pickled = pickle.load(open(filepath, 'rb'), encoding='latin1')
    else:
        pickled = pickle.load(open(filepath, 'rb'), encoding='utf-8')
    return pickled


def pad_or_truncate(some_list, target_len):
    """Pad or truncate a list upto given target length

    Parameters
    ----------
    some_list : list
                Input list
    target_length : int
                    Final length of list

    If being extended, returns list padded with NAs.
    """
    return some_list[:target_len] + [np.nan] * (target_len - len(some_list))


def pad_five_prime_or_truncate(some_list, offset_5p, target_len):
    """Pad first the 5prime end and then the 3prime end or truncate

    Parameters
    ----------
    some_list : list
                Input list
    offset_5p : int
                5' offset
    target_length : int
                    Final length of list

    If being extended, returns list padded with NAs.
    """
    some_list = list(some_list)
    padded_5p = [np.nan] * offset_5p + some_list
    return padded_5p[:target_len] + [np.nan] * (target_len - len(padded_5p))


def codon_to_anticodon(codon):
    """Codon to anticodon.

    Parameters
    ----------
    codon : string
            Input codon
    """
    pairs = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N'}
    return ''.join(pairs[c] for c in codon)[::-1]


def collapse_bed_intervals(intervals,
                           chromosome_lengths=None,
                           offset_5p=0,
                           offset_3p=0):
    """Collapse intervals into non overlapping manner


    # NOTE
    # TODO : This function has a subtle bug that it will be offset by 1
    # position when the gene is on negative strand
    # So essentially if you have CDS on a negative strand
    # The first position should be discarded
    # Similary for the last position in the gene on + strand
    # you have an extra position in the end

    Parameters
    ----------
    intervals : list of tuples
                Like [('chr1', 310, 320, '+'), ('chr1', 321, 330, '+')]
    chromosome_lengths : dict
                         A map of each chromosome'e length
                         Only used with offset_3p, offset_5p>0
    offset_5p : int (positive)
                Number of bases to count upstream (5')
    offset_3p : int (positive)
                Number of bases to count downstream (3')

    Returns
    -------
    interval_combined : list of tuples
                        A collapsed version of interval
                        This is useful when the annotations are overlapping.
                        Example:
                        chr1 310 320 gene1 +
                        chr1 319 324 gene1 +
                        Returns:
                        chr1 310 324 gene1 +

    intervals_for_fasta_read : list of tuples
                               This list can be used to directly fetch
                               fasta from pyfaidx.
                               NOTE: DO NOT do offset adjustments
                               as they are already adjusted for pyfaidx format
                               (1-end both start and end)

    gene_offset_5p, gene_offset_3 : in
                                    Gene wise offsets.
                                    This might be different from offset_5p in cases where
                                    `offset_5p` leads to a negative coordinate
    """
    chrom = intervals[0][0]
    strand = intervals[0][3]
    chroms = list(set([i[0] for i in intervals]))
    strands = list(set([i[3] for i in intervals]))

    if len(chroms) != 1:
        sys.stderr.write('Error chromosomes should be unique')
        return
    if len(strands) != 1:
        sys.stderr.write('Error strands should be unique')
        return

    intervals_list = [list(element) for element in list(intervals)[:]]

    first_interval = intervals_list[0]
    last_interval = intervals_list[-1]
    if offset_5p != 0 and offset_3p != 0:
        chrom_length = chromosome_lengths[str(first_interval[0])]
    else:
        chrom_length = np.inf
    # Need to convert to list instead frm tuples
    # TODO fix this?
    # intervals = list(map(list, list(intervals)))
    if strand == '+':
        # For positive strand shift
        # start codon position first_interval[1] by -offset
        if first_interval[1] - offset_5p >= 0:
            first_interval[1] = first_interval[1] - offset_5p
            gene_offset_5p = offset_5p
        else:
            sys.stderr.write('Cannot offset beyond 0 for interval: {}. \
                Set to start of chromsome.\n'.format(first_interval))
            # Reset offset to minimum possible
            gene_offset_5p = first_interval[1]
            first_interval[1] = 0

        if (last_interval[2] + offset_3p <= chrom_length):
            last_interval[2] = last_interval[2] + offset_3p
            gene_offset_3p = offset_3p
        else:
            sys.stderr.write('Cannot offset beyond 0 for interval: {}. \
                             Set to end of chromsome.\n'.format(last_interval))
            gene_offset_3p = chrom_length - last_interval[2]
            # 1-end so chrom_length
            last_interval[2] = chrom_length
    else:
        # Else shift cooridnate of last element in intervals stop by + offset
        if (last_interval[2] + offset_5p <= chrom_length):
            last_interval[2] = last_interval[2] + offset_5p
            gene_offset_5p = offset_5p
        else:
            sys.stderr.write('Cannot offset beyond 0 for interval: {}. \
                             Set to end of chromsome.\n'.format(last_interval))
            gene_offset_5p = chrom_length - last_interval[2]
            # 1-end so chrom_length
            last_interval[2] = chrom_length
        if first_interval[1] - offset_3p >= 0:
            first_interval[1] = first_interval[1] - offset_3p
            gene_offset_3p = offset_3p
        else:
            sys.stderr.write('Cannot offset beyond 0 for interval: {}. \
                Set to start of chromsome.\n'.format(first_interval))
            # Reset offset to minimum possible
            gene_offset_5p = first_interval[1]
            first_interval[1] = 0

    intervals = [tuple(element) for element in intervals_list]

    interval_coverage_list = []
    for index, interval in enumerate(intervals):
        strand = interval[3]
        if strand == '+':
            series_range = list(range(interval[1], interval[2]))
        elif strand == '-':
            series_range = list(range(interval[2], interval[1], -1))

        series = pd.Series(series_range, index=series_range)
        interval_coverage_list.append(series)

    if len(interval_coverage_list) == 0:
        # Some genes might not be present in the bigwig at all
        sys.stderr.write('Got empty list! intervals  for chr : {}\n'.format(
            first_interval[0]))
        return (pd.Series([]), pd.Series([]), pd.Series([]), 0, 0)

    interval_combined = interval_coverage_list[0]
    for interval in interval_coverage_list[1:]:
        interval_combined = interval_combined.combine_first(interval)
    interval_index = np.arange(len(interval_combined)) - gene_offset_5p
    index_to_genomic_pos_map = pd.Series(
        interval_combined.index.tolist(), index=interval_index)
    """
    intervals_for_fasta_read = []
    for pos in index_to_genomic_pos_map.values:
        # we use 1-based indexing (both start and end) for fetching fasta
        intervals_for_fasta_read.append((chrom, pos + 1, pos + 1, strand))
    """
    interval_combined = interval_combined.reset_index(drop=True)
    interval_combined = interval_combined.rename(lambda x: x - gene_offset_5p)
    interval_zero_start = list_to_ranges(
        interval_combined.sort_values().astype(int).values.tolist())
    interval_one_start = list_to_ranges(
        (1 + interval_combined).sort_values().astype(int).values.tolist())
    query_intervals = []
    for index, (start, stop) in enumerate(interval_zero_start):
        # our intervals were already 0-based start and 1-based end
        # so we want to retain that
        # it is tricky
        # but the start and end were both transformed to 0-based since we used range before
        query_intervals.append((chrom, start, stop + 1, strand))
    fasta_onebased_intervals = []
    for index, (start, stop) in enumerate(interval_one_start):
        # our intervals were already 0-based start and 1-based end
        # so we want to retain that
        # it is tricky
        # but the start and end were both transformed to 0-based since we used range before
        fasta_onebased_intervals.append((chrom, start, stop, strand))
    return query_intervals, fasta_onebased_intervals, gene_offset_5p, gene_offset_3p
