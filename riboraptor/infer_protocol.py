from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pysam
from .count import _is_read_uniq_mapping, _create_bam_index
from .helpers import read_refseq_bed
from collections import Counter
import numpy as np


def infer_protocol(bam, refseq, n_reads=20000):
    """Infer strandedness protocol given a bam file

    Parameters
    ----------
    bam: string
         Path to bam file
    refseq: string
                Path to refseq bed file
    n_reads: int
             Number of reads to use (downsampled)


    Returns
    -------
    protocol: string
              unstranded/forward/reverse
    forward_mapped_reads: float
          Proportion of reads of type + mapping to + (++) or - mapping to - (--)
    reverse_mapped_reads: float
          Proportion of reads of type + mapping to - (+-) or - mapping to + (-+)
    """
    iteration = 0
    bam = pysam.AlignmentFile(bam)
    refseq = read_refseq_bed(refseq)
    strandedness = Counter()
    for read in bam:
        if not _is_read_uniq_mapping(read):
            continue
        if read.is_reverse:
            mapped_strand = '-'
        else:
            mapped_strand = '+'
        mapped_start = read.reference_start
        mapped_end = read.reference_end
        chrom = read.reference_name
        gene_strand = list(set(refseq[chrom].find(mapped_start, mapped_end)))
        if len(gene_strand) != 1:
            # Filter out genes with ambiguous strand info
            # (those) that have a tx_start on opposite strands
            continue
        gene_strand = gene_strand[0]
        strandedness['{}{}'.format(mapped_strand, gene_strand)] += 1
        iteration += 1
        if iteration >= n_reads:
            break
    ## Add pseudocounts
    strandedness['++'] += 1
    strandedness['--'] += 1
    strandedness['+-'] += 1
    strandedness['-+'] += 1

    total = sum(strandedness.values())
    forward_mapped_reads = (strandedness['++'] + strandedness['--']) / total
    reverse_mapped_reads = (strandedness['-+'] + strandedness['+-']) / total
    ratio = forward_mapped_reads / reverse_mapped_reads
    if np.isclose([ratio], [1]):
        return 'unstranded', forward_mapped_reads, reverse_mapped_reads
    elif forward_mapped_reads >= 0.5:
        return 'forward', forward_mapped_reads, reverse_mapped_reads
    else:
        return 'reverse', forward_mapped_reads, reverse_mapped_reads
