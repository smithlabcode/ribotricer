"""Plotting methods."""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

from collections import Counter
from itertools import cycle
from itertools import islice

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FormatStrFormatter

import numpy as np
import pandas as pd
import seaborn as sns
import pycwt as wavelet

from .utils.helpers import identify_peaks
from .utils.helpers import millify
from .utils.helpers import round_to_nearest

__FRAME_COLORS__ = ['#1b9e77', '#d95f02', '#7570b3']


def setup_plot():
    """Setup plotting defaults"""
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


def setup_axis(ax, axis='x', majorticks=5,
               minorticks=1, xrotation=0, yrotation=0):
    """Setup axes defaults

    Parameters
    ----------

    ax : matplotlib.Axes

    axis : str
           Setup 'x' or 'y' axis
    majorticks : int
                 Length of interval between two major ticks
    minorticks : int
                 Length of interval between two major ticks
    xrotation : int
                Rotate x axis labels by xrotation degrees
    yrotation : int
                Rotate x axis labels by xrotation degrees
    """
    major_locator = MultipleLocator(majorticks)
    major_formatter = FormatStrFormatter('%d')
    minor_locator = MultipleLocator(minorticks)
    if axis == 'x':
        ax.xaxis.set_major_locator(major_locator)
        ax.xaxis.set_major_formatter(major_formatter)
        ax.xaxis.set_minor_locator(minor_locator)
    elif axis == 'y':
        ax.yaxis.set_major_locator(major_locator)
        ax.yaxis.set_major_formatter(major_formatter)
        ax.yaxis.set_minor_locator(minor_locator)
    elif axis == 'both':
        setup_axis(ax, 'x', majorticks, minorticks, xrotation, yrotation)
        setup_axis(ax, 'y', majorticks, minorticks, xrotation, yrotation)
    ax.tick_params(which='major', width=2, length=10)
    ax.tick_params(which='minor', width=1, length=6)


def plot_fragment_dist(fragment_lengths, ax=None, millify_labels=True, **kwargs):
    """Plot fragment length distribution.

    Parameters
    ----------
    fragment_lengths : array_like
                     Array of fragment lengths

    ax : matplotlib.Axes
        Axis object
    millify_labels : bool
                     True if labels should be formatted to read millions/trillions etc

    """
    # setup_plot()
    fig = None
    if ax is None:
        fig, ax = plt.subplots()
    setup_axis(ax, **kwargs)
    if isinstance(fragment_lengths, Counter):
        fragment_lengths = pd.Series(fragment_lengths)
        fragment_lengths_counts = fragment_lengths.sort_index()
    else:
        fragment_lengths = pd.Series(fragment_lengths)
        fragment_lengths_counts = fragment_lengths.value_counts().sort_index()

    ax.bar(fragment_lengths_counts.index, fragment_lengths_counts)
    ax.set_xlim(min(fragment_lengths_counts.index) - 0.5,
                round_to_nearest(max(fragment_lengths_counts.index), 10) + 0.5)
    if millify_labels:
        ax.set_yticklabels(list(map(lambda x: millify(x), ax.get_yticks())))
    sns.despine(trim=True, offset=20)
    return ax, fig


def plot_framewise_dist(counts, fragment_len_range, ax=None):
    """Plot framewise distribution of fragments.

    Parameters
    ----------
    counts : Series
            A series with position as index and value as counts

    fragment_len_range: int or range
        Range of fragment lengths to average counts over

    ax : matplotlib.Axes
        Default none

    """
    setup_plot()
    assert isinstance(counts, pd.Series)
    if ax is None:
        _, ax = plt.subplots()
    setup_axis(ax)
    ax.set_ylabel('Number of reads')
    ax.set_xlim(min(counts.index) - 0.6,
                round_to_nearest(max(counts.index), 10) + 0.6)
    barlist = ax.bar(counts.index, counts.values)
    barplot_colors = list(
        islice(cycle(__FRAME_COLORS__), None, len(counts.index)))
    for index, cbar in enumerate(barlist):
        cbar.set_color(barplot_colors[index])
    ax.legend((barlist[0], barlist[1], barlist[2]),
              ('Frame 1', 'Frame 2', 'Frame 3'))

    sns.despine(trim=True, offset=20)
    return ax


def plot_read_counts(counts, ax=None, marker=False, color='royalblue',
                     label=None, millify_labels=True, identify_peak=True, **kwargs):
    """Plot RPF density around start/stop codons.

    Parameters
    ----------
    counts : Series/Counter
             A series with coordinates as index and counts as values
    ax : matplotlib.Axes
         Axis to create object on
    marker : string
             'o'/'x'
    color : string
            Line color
    label : string
            Label (useful only if plotting multiple objects on same axes)
    millify_labels : bool
                     True if labels should be formatted to read millions/trillions etc

    """
    # setup_plot()
    if isinstance(counts, Counter):
        counts = pd.Series(counts)
    fig = None
    if ax is None:
        fig, ax = plt.subplots()
    setup_axis(ax, **kwargs)
    ax.set_ylabel('Number of reads')
    if not marker:
        ax.plot(counts.index, counts.values, color=color,
                linewidth=2, label=label)
    else:
        ax.plot(counts.index, counts.values, color=color,
                marker='o', linewidth=2, label=label)
    # ax.set_xlim(round_to_nearest(ax.get_xlim()[0], 50) - 0.6,
    #            round_to_nearest(ax.get_xlim()[1], 50) + 0.6)
    peak = None
    if identify_peak:
        peak = identify_peaks(counts)
        ax.axvline(x=peak, color='r', linestyle='dashed')
        ax.text(peak + 0.5, ax.get_ylim()[1]
                * 0.9, '{}'.format(peak), color='r')

    if millify_labels:
        ax.set_yticklabels(list(map(lambda x: millify(x), ax.get_yticks())))
    ax.set_xlim(min(counts.index) - 0.5,
                round_to_nearest(max(counts.index), 10) + 0.5)
    sns.despine(trim=True, offset=10)
    return ax, fig, peak


def create_wavelet(data, ax):
    t = data.index

    N = len(data.index)
    p = np.polyfit(data.index, data, 1)
    data_notrend = data - np.polyval(p, data.index)
    std = data_notrend.std()  # Standard deviation
    var = std**2  # Variance
    data_normalized = data_notrend / std  # Normalized dataset

    mother = wavelet.Morlet(6)
    dt = 1
    s0 = 2 * dt  # Starting scale, in this case 2 * 0.25 years = 6 months
    dj = 1 / 12  # Twelve sub-octaves per octaves
    J = 7 / dj  # Seven powers of two with dj sub-octaves
    alpha, _, _ = wavelet.ar1(data)  # Lag-1 autocorrelation for red noise

    wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(
        data_normalized, dt=dt, dj=dj, s0=s0, J=J, wavelet=mother)
    iwave = wavelet.icwt(wave, scales, dt, dj, mother) * std

    power = (np.abs(wave))**2
    fft_power = np.abs(fft)**2
    period = 1 / freqs

    power /= scales[:, None]
    signif, fft_theor = wavelet.significance(
        1.0, dt, scales, 0, alpha, significance_level=0.95, wavelet=mother)
    sig95 = np.ones([1, N]) * signif[:, None]
    sig95 = power / sig95

    glbl_power = power.mean(axis=1)
    dof = N - scales  # Correction for padding at edges
    glbl_signif, tmp = wavelet.significance(
        var,
        dt,
        scales,
        1,
        alpha,
        significance_level=0.95,
        dof=dof,
        wavelet=mother)

    levels = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16]
    ax.contourf(
        t,
        np.log2(period),
        np.log2(power),
        np.log2(levels),
        extend='both',
        cmap=plt.cm.viridis)
    extent = [t.min(), t.max(), 0, max(period)]
    ax.contour(
        t,
        np.log2(period),
        sig95, [-99, 1],
        colors='k',
        linewidths=2,
        extent=extent)
    ax.fill(
        np.concatenate([t, t[-1:] + dt, t[-1:] + dt, t[:1] - dt, t[:1] - dt]),
        np.concatenate([
            np.log2(coi), [1e-9],
            np.log2(period[-1:]),
            np.log2(period[-1:]), [1e-9]
        ]),
        'k',
        alpha=0.3,
        hatch='x')
    ax.set_title('Wavelet Power Spectrum')
    ax.set_ylabel('Frequency')
    Yticks = 2**np.arange(0, np.ceil(np.log2(period.max())))
    ax.set_yticks(np.log2(Yticks))
    ax.set_yticklabels(np.round(1 / Yticks, 3))

    return (iwave, period, power, sig95, coi)
