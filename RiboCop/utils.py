"""Utilities for analysis
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import warnings

from collections import Counter
from collections import defaultdict

def parse_ccds(annotation, orfs, saveto):
    """
    annotation: str
                Path for annotation files of putative ORFs
    orfs: str
          Path for translating ORFs
    saveto: str
          output file name
    """
    ccds = defaultdict(list)
    with open(annotation, 'r') as anno:
        for line in anno:
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

    orfs = {}
    with open(orfs, 'r') as orf:
        for line in orf:
            if not line:
                print('orf line cannot be empty')
                return None
            fields = line.split('\t')
            if len(fields) != 5:
                print('unexpected number of columns found for orf file')
                return None
            oid = fields[0]
            count = fields[2]
            corr = fields[3]
            pval = fields[4]
            orfs[oid] = (count, corr, pval)

    to_write = 'Gene_ID\tCount\tPeriodicity\tPval\n'
    for gid in ccds:
        count, corr, pval = (0, 0, 0)
        for oid in ccds[gid]:
            t_cnt, t_corr, t_pval = orfs[oid]
            if t_corr >= corr:
                count, corr, pval = (t_cnt, t_corr, t_pval)
        to_write += '{}\t{}\t{}\t{}\n'.format(gid, count, corr, pval)

    with open(saveto, 'w') as output:
        output.write(to_write)


