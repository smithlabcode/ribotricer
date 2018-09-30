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
                    print('unexpected number of columns found for annotation file')
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


