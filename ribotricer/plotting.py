"""Plotting functions."""
# Part of ribotricer software
#
# Copyright (C) 2019 Wenzheng Li, Saket Choudhary and Andrew D Smith
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.use("Agg")


def plot_read_lengths(read_lengths, prefix):
    """
    Parameters
    ----------
    read_lengths: dict
                  key is the length, value is the number of reads
    prefix: str
            prefix for the output file
    """
    # print('plotting read length distribution...')
    fig, ax = plt.subplots()
    x = sorted(read_lengths.keys())
    y = [read_lengths[i] for i in x]
    ax.bar(x, y)
    ax.set_xlabel("Read length")
    ax.set_ylabel("Number of reads")
    ax.set_title("Read length distribution")
    fig.tight_layout()
    fig.savefig("{}_read_length_dist.pdf".format(prefix))
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
    # print('plotting metagene profiles...')
    total_reads = sum(read_lengths.values())
    frame_colors = ["#fc8d62", "#66c2a5", "#8da0cb"]
    with PdfPages("{}_metagene_plots.pdf".format(prefix)) as pdf:
        for length in sorted(metagenes):
            # TODO: This only consider the 5' end, should be generalized to 3'
            metagene_cov_start, metagene_cov_stop, coh, valid, _, _ = metagenes[length]
            if len(metagene_cov_start) != 0:
                min_index = min(metagene_cov_start.index.tolist())
                max_index = max(metagene_cov_start.index.tolist())
                start_offset = min(offset, max_index)
                metagene_cov_start = metagene_cov_start[
                    np.arange(min_index, start_offset)
                ]
                x = np.arange(min_index, start_offset)
                colors = np.tile(frame_colors, len(x) // 3 + 1)
                xticks = np.arange(min_index, start_offset, 20)
                ratio = read_lengths[length] / total_reads
                fig, (ax, ax2) = plt.subplots(nrows=2, ncols=1)
                ax.vlines(
                    x, ymin=np.zeros(len(x)), ymax=metagene_cov_start, colors=colors
                )
                ax.tick_params(axis="x", which="both", top=False, direction="out")
                ax.set_xticks(xticks)
                ax.set_xlim((min_index, start_offset))
                ax.set_xlabel("Distance from start codon (nt)")
                ax.set_ylabel("Normalized mean reads")
                ax.set_title(
                    ("{} nt reads, proportion: {:.2%}\nphase_score: {:.2}").format(
                        length, ratio, coh
                    )
                )

                # plot distance from stop codon
                min_index = min(metagene_cov_stop.index.tolist())
                max_index = max(metagene_cov_stop.index.tolist())
                stop_offset = max(-offset, min_index)
                metagene_cov_stop = metagene_cov_stop[np.arange(stop_offset, max_index)]
                x = np.arange(stop_offset, max_index)
                colors = np.tile(frame_colors, len(x) // 3 + 1)
                xticks = np.arange(stop_offset, max_index, 20)
                ax2.vlines(
                    x, ymin=np.zeros(len(x)), ymax=metagene_cov_stop, colors=colors
                )
                ax2.tick_params(axis="x", which="both", top=False, direction="out")
                ax2.set_xticks(xticks)
                ax2.set_xlim((stop_offset, max_index))
                ax2.set_xlabel("Distance from stop codon (nt)")
                ax2.set_ylabel("Normalized mean reads")

                fig.tight_layout()
                pdf.savefig(fig)
                plt.close()
