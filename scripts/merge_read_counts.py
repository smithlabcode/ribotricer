"""Utilities for creating read count tables
"""

from os import listdir
from os.path import isfile, join
import sys
import pandas as pd
import numpy as np
import collections
from functools import reduce

def main():
    project = sys.argv[1]
    fname = sys.argv[2]
    with open(fname) as f:
        content = f.readlines()
    samples = [x.strip() for x in content]

    regions = ['UTR5', 'CDS', 'UTR3']
    read_counts = []
    root = '/staging/as/wenzhenl/re-ribo-analysis/{}/mapped/read_counts/{}'

    for region in regions:
        path = root.format(project, region)
        region_tables = []
        for sample in samples:
            df = pd.read_table(join(path, sample + '_read_counts.tsv'))
            df = df[['gene_name', 'count']]
            df.rename(columns = {'count': sample})
            region_tables.append(df)
        df_region = reduce(lambda left, right: 
                           pd.merge(left, right,
                                    on='gene_name'),
                           region_tables)
        saveto = join(root.format(project, region), region + "cnt_table.csv")
        df_region.to_csv(saveto, sep=str('\t'), index=False)

if __name__ == "__main__":
    main()

