"""Plotting functions."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from .statistics import coherence


def plot_read_lengths(read_lengths, prefix):
    """
    Parameters
    ----------
    read_lengths: dict
                  key is the length, value is the number of reads
    prefix: str
            prefix for the output file
    """
    print('plotting read length distribution...')
    fig, ax = plt.subplots()
    x = sorted(read_lengths.keys())
    y = [read_lengths[i] for i in x]
    ax.bar(x, y)
    ax.set_xlabel('Read length')
    ax.set_ylabel('Number of reads')
    ax.set_title('Read length distribution')
    fig.tight_layout()
    fig.savefig('{}_read_length_dist.pdf'.format(prefix))
    plt.close()


def plot_metagene(metagenes, read_lengths, prefix, offset=200):
    """
    Parameters
    ----------
    metagenes: dict
               key is the length, value is the metagene coverage
    read_lengths: dict
                  key is the length, value is the number of reads
    prefix: str
            prefix for the output file
    """
    print('plotting metagene profiles...')
    total_reads = sum(read_lengths.values())
    with PdfPages('{}_metagene_plots.pdf'.format(prefix)) as pdf:
        for length in sorted(metagenes):
            metagene_cov_start, metagene_cov_stop = metagenes[length]
            if len(metagene_cov_start) == 0:
                continue
            corr, pval, nonzero = coherence(metagene_cov_start.values)
            min_index = min(metagene_cov_start.index.tolist())
            max_index = max(metagene_cov_start.index.tolist())
            offset = min(offset, max_index)
            metagene_cov_start = metagene_cov_start[np.arange(
                min_index, offset)]
            x = np.arange(min_index, offset)
            colors = np.tile(['r', 'g', 'b'], len(x) // 3 + 1)
            xticks = np.arange(min_index, offset, 20)
            ratio = read_lengths[length] / total_reads
            fig, (ax, ax2) = plt.subplots(nrows=2, ncols=1)
            ax.vlines(
                x,
                ymin=np.zeros(len(x)),
                ymax=metagene_cov_start,
                colors=colors)
            ax.tick_params(axis='x', which='both', top='off', direction='out')
            ax.set_xticks(xticks)
            ax.set_xlim((min_index, offset))
            ax.set_xlabel('Distance from start codon (nt)')
            ax.set_ylabel('Normalized mean reads')
            ax.set_title((
                '{} nt reads, proportion: {:.2%}\nPeriodicity: {:.2}, pval: {:.6}'
            ).format(length, ratio, corr, pval))

            ### plot distance from stop codon
            min_index = min(metagene_cov_stop.index.tolist())
            max_index = max(metagene_cov_stop.index.tolist())
            offset = max(-offset, min_index)
            metagene_cov_stop = metagene_cov_stop[np.arange(offset, max_index)]
            x = np.arange(offset, max_index)
            colors = np.tile(['r', 'g', 'b'], len(x) // 3 + 1)
            xticks = np.arange(offset, max_index, 20)
            ax2.vlines(
                x,
                ymin=np.zeros(len(x)),
                ymax=metagene_cov_stop,
                colors=colors)
            ax2.tick_params(axis='x', which='both', top='off', direction='out')
            ax2.set_xticks(xticks)
            ax2.set_xlim((offset, max_index))
            ax2.set_xlabel('Distance from stop codon (nt)')
            ax2.set_ylabel('Normalized mean reads')

            fig.tight_layout()
            pdf.savefig(fig)
            plt.close()
