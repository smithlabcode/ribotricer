"""Utilities for analysis
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import warnings

from collections import Counter
from collections import defaultdict
from tqdm import *
from .common import cal_periodicity
from .common import coherence


def parse_ccds(annotation, orfs, saveto):
    """
    annotation: str
                Path for annotation files of putative ORFs
    orfs: str
          Path for translating ORFs
    saveto: str
          output file name
    """
    anno_oids = []
    real_oids = []
    ccds = defaultdict(list)
    with open(annotation, 'r') as anno:
        total_lines = len(['' for line in anno])
    with open(annotation, 'r') as anno:
        with tqdm(total=total_lines) as pbar:
            header = True
            for line in anno:
                pbar.update()
                if header:
                    header = False
                    continue
                if not line:
                    print('annotation line cannot be empty')
                    return None
                fields = line.split('\t')
                if len(fields) != 13:
                    print(
                        'unexpected number of columns found for annotation file'
                    )
                    return None
                oid = fields[0]
                gid = fields[4]
                ccds[gid].append(oid)
                anno_oids.append(oid)

    ccds_orfs = {}
    with open(orfs, 'r') as orf:
        total_lines = len(['' for line in orf])
    with open(orfs, 'r') as orf:
        with tqdm(total=total_lines) as pbar:
            header = True
            for line in orf:
                pbar.update()
                if header:
                    header = False
                    continue
                if not line:
                    print('orf line cannot be empty')
                    return None
                fields = line.split('\t')
                if len(fields) != 5:
                    print('unexpected number of columns found for orf file')
                    return None
                oid = fields[0]
                count = int(fields[2])
                corr = float(fields[3])
                pval = float(fields[4])
                ccds_orfs[oid] = (count, corr, pval)
                real_oids.append(oid)

    rename = {x: y for (x, y) in zip(anno_oids, real_oids)}
    to_write = 'Gene_ID\tCount\tPeriodicity\tPval\n'
    n_genes = 0
    for gid in ccds:
        n_genes += 1
        count, corr, pval = (0, 0, 1)
        for oid in ccds[gid]:
            oid = rename[oid]
            t_cnt, t_corr, t_pval = ccds_orfs[oid]
            if t_corr >= corr:
                count, corr, pval = (t_cnt, t_corr, t_pval)
        to_write += '{}\t{}\t{}\t{}\n'.format(gid, count, corr, pval)
    print(n_genes)

    with open(saveto, 'w') as output:
        output.write(to_write)

def test_periodicity(orf_file, prefix, method):

    print('testing method: {}'.format(method))
    print('exporting coverages for all ORFs...')
    to_write = 'ORF_ID\tcoverage\tcount\tlength\tnonzero\tperiodicity\tpval\n'
    with open(orf_file, 'r') as orf:
        total_lines = len(['' for line in orf])
    with open(orf_file, 'r') as  orf:
        header = True
        with tqdm(total=total_lines) as pbar:
            for line in orf:
                pbar.update()
                if header:
                    header = False
                    continue
                fields = line.split('\t')
                oid = fields[0]
                cov = fields[1]
                cov = cov[1:-1]
                cov = [int(x) for x in cov.split(', ')]
                count = sum(cov)
                length = len(cov)
                if len(cov) < 60:
                    corr, pval, nonzero = (0, 1, 0)
                else:
                    corr, pval, nonzero = cal_periodicity(cov)

                to_write += '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    oid, cov, count, length, nonzero, corr, pval)
    with open('{}_translating_ORFs_{}.tsv'.format(prefix, method), 'w') as output:
        output.write(to_write)

def benchmark(rna_file, ribo_file, prefix, cutoff=5):

    rna = {}
    ribo = {}

    print('reading RNA profiles')
    with open(rna_file, 'r') as orf:
        total_lines = len(['' for line in orf])
    with open(rna_file, 'r') as  orf:
        with tqdm(total=total_lines) as pbar:
            for line in orf:
                pbar.update()
                chrom, start, end, cat, gid, strand, cov = line.strip().split('\t')
                cov = [int(x) for x in cov.strip().split()]
                if strand == '-':
                    cov.reverse()
                ID = '_'.join([chrom, start, end, cat, gid])
                rna[ID] = cov

    print('reading Ribo profiles')
    with open(ribo_file, 'r') as orf:
        total_lines = len(['' for line in orf])
    with open(ribo_file, 'r') as  orf:
        with tqdm(total=total_lines) as pbar:
            for line in orf:
                pbar.update()
                chrom, start, end, cat, gid, strand, cov = line.strip().split('\t')
                cov = [int(x) for x in cov.strip().split()]
                if strand == '-':
                    cov.reverse()
                ID = '_'.join([chrom, start, end, cat, gid])
                ribo[ID] = cov

    to_write = 'ID\tmy_ribo\tmy_rna\tribo_coh\trna_coh\tribo_cov\trna_cov\n'
    common_ids = set(ribo.keys()) & set(rna.keys())
    for ID in tqdm(common_ids):
        if sum(rna[ID]) < cutoff or sum(ribo[ID]) < cutoff:
            continue
        if len(ribo[ID]) < 10:
            continue
        rna_coh, rna_pval, rna_valid = coherence(rna[ID])
        rna_cov = rna_valid / len(rna[ID])
        if rna_valid // 3 < 4:
            continue
        ribo_coh, ribo_pval, ribo_valid = coherence(ribo[ID])
        ribo_cov = ribo_valid / len(ribo[ID])
        if ribo_valid // 3 < 4:
            continue

        to_write += '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    ID, ribo_pval, rna_pval, ribo_coh, rna_coh, ribo_cov, rna_cov)
    with open('{}_results.txt'.format(prefix), 'w') as output:
        output.write(to_write)
