"""Plotting methods."""
from itertools import cycle
from itertools import islice
import matplotlib.pyplot as plt

import seaborn as sns
import pandas as pd
from matplotlib.ticker import AutoMinorLocator

from .helpers import round_to_nearest

from .count_utilities import summarize_counts

def setup_plot():
    plt.rcParams['savefig.dpi'] = 120
    plt.rcParams['figure.dpi'] = 120
    plt.rcParams['figure.autolayout'] = False
    plt.rcParams['figure.figsize'] = 12, 8
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['axes.titlesize'] = 20
    plt.rcParams['font.size'] = 10
    plt.rcParams['lines.linewidth'] = 2.0
    plt.rcParams['lines.markersize'] = 8
    plt.rcParams['legend.fontsize'] = 14

    sns.set_style('white')
    sns.set_context('paper', font_scale=2)

def setup_axis(ax, minorticks=10):

    minorLocator = AutoMinorLocator(minorticks)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.tick_params(which='major', width=2)
    ax.tick_params(which='minor', width=1)
    ax.tick_params(which='major', length=10)
    ax.tick_params(which='minor', length=6)
#__FRAME_COLORS__ = ['#2ecc71', '#3498db', '#e74c3c']
__FRAME_COLORS__ = ['#1b9e77',  '#d95f02', '#7570b3']


def plot_framewise_dist(all_frag_counts, fragment_len_range, ax=None):
    """Plot framewise distribution of fragments.

    Parameters
    ----------
    all_frag_counts : dict or Series
            A pandas Seires object with the position as key and each values as a dict of fragment length and read counts

    fragment_len_range: int or range
        Range of fragment lengths to average counts over

    ax : matplotlib.Axes
        Default none

    """
    setup_plot()
    counts = summarize_counts(all_frag_counts, fragment_len_range)


    assert isinstance(counts, pd.Series)

    if ax is None:
        fig, ax = plt.subplots()
    setup_axis(ax)
    ax.set_ylabel('Number of reads')
    ax.set_xlim(min(counts.index) - 0.6,
                round_to_nearest(max(counts.index), 10) + 0.6)
    barlist = ax.bar(counts.index, counts.values)
    barplot_colors = list(islice(cycle(__FRAME_COLORS__), None, len(counts.index)))
    for index, bar in enumerate(barlist):
        bar.set_color(barplot_colors[index])
    ax.legend((barlist[0], barlist[1], barlist[2]),
              ('Frame 1', 'Frame 2', 'Frame 3'))

    sns.despine(trim=True, offset=20)
    return ax


def plot_fragment_dist(fragment_lengths, ax=None):
    """Plot fragment length distribution.

    Parameters
    ----------
    fragment_lengths : array_like
                     Array of fragment lengths

    ax : matplotlib.Axes
        Axis object

    """
    setup_plot()
    if ax is None:
        fig, ax = plt.subplots()
    setup_axis(ax, 5)
    fragment_lengths = pd.Series(fragment_lengths)
    fragment_lengths_counts = fragment_lengths.value_counts().sort_index()

    ax.bar(fragment_lengths_counts.index, fragment_lengths_counts)
    ax.set_xlim(min(fragment_lengths_counts.index) - 0.5,
                round_to_nearest(max(fragment_lengths_counts.index), 10) + 0.5)
    sns.despine(trim=True, offset=20)
    return ax


def plot_continuous(all_frag_counts, fragment_len_range, ax=None, half_window=30):
    setup_plot()
    counts = summarize_counts(all_frag_counts, fragment_len_range)

    counts = counts[range(-half_window, half_window)]
    assert isinstance(counts, pd.Series)
    setup_plot()
    if ax is None:
        fig, ax = plt.subplots()

    setup_axis(ax)

    ax.set_ylabel('Number of reads')
    ax.set_xlim(min(counts.index) - 0.6, round_to_nearest(max(counts.index), 10) + 0.6)
    ax.plot(counts.index, counts.values, color='royalblue', marker='o', linewidth=2)#'#3498db')
    sns.despine(trim=True, offset=20)
    return ax
