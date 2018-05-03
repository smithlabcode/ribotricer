from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
import glob
import re
import pickle
import os

import pandas as pd


def summary_starlogs_over_runs(directory, list_of_srr):
    df = pd.DataFrame()
    files_not_found = []
    for run in list_of_srr:
        if not os.path.isfile(os.path.join(directory, run + '.tsv')):
            files_not_found.append(run)
            continue
        temp_df = pd.read_table(os.path.join(directory, run + '.tsv'))
        df = pd.concat([df, temp_df])
    return df, files_not_found


def load_tpm(path):
    df = pd.read_table(path, names=['gene_id', 'tpm']).set_index('gene_id')
    return df


def get_cell_line_or_tissue(row):
    if str(row['cell_line']).strip() and str(row['cell_line']).strip() != 'nan':
        return '{}-{}-{}'.format(row['cell_line'], row['study_accession'],
                                 row['experiment_accession'])
    if str(row['tissue']).strip() and str(row['tissue']).strip() != 'nan':
        return '{}-{}-{}'.format(row['tissue'], row['study_accession'],
                                 row['experiment_accession'])
    if str(row['source_name']).strip(
    ) and str(row['source_name']).strip() != 'nan':
        return '{}-{}-{}'.format(row['source_name'], row['study_accession'],
                                 row['experiment_accession'])
    if row['study_accession'].strip() == 'SRP052229':
        print(row)
    return '{}-{}-{}'.format(row['source_name'], row['study_accession'],
                             row['experiment_accession'])


def determine_cell_type(sample_attribute):
    sample_attribute = str(sample_attribute)
    if 'cell line:' in sample_attribute:
        x = re.search(r'cell line: \w+', sample_attribute)
        return x.group(0).strip('cell line: ').rstrip(' ').upper()
    if 'cell_line:' in sample_attribute:
        x = re.search(r'cell_line: \w+', sample_attribute)
        return x.group(0).strip('cell_line: ').rstrip(' ').upper()
    if 'cell-line:' in sample_attribute:
        x = re.search(r'cell-line: \w+', sample_attribute)
        return x.group(0).strip('cell-line: ').rstrip(' ').upper()
    if 'cell_type:' in sample_attribute:
        x = re.search(r'cell_type: \w+', sample_attribute)
        return x.group(0).strip('cell_type: ').rstrip(' ').upper()
    if 'source_name:' in sample_attribute:
        x = re.search(r'source_name: \w+', sample_attribute)
        return x.group(0).strip('source_name: ').rstrip(' ').upper()
    else:
        #pass
        print('Found {}'.format(sample_attribute))
        return np.nan


def get_tissue_type(sample_attribute):
    sample_attribute = str(sample_attribute)
    if 'tissue: ' in sample_attribute:
        x = re.search(r'tissue: \w+', sample_attribute)
        return x.group(0).strip('tissue: ').rstrip(' ').lower()
    else:
        print('Found {}'.format(sample_attribute))
        return np.nan


def get_strain_type(sample_attribute):
    sample_attribute = str(sample_attribute)
    if 'strain: ' in sample_attribute:
        x = re.search(r'strain: \w+', sample_attribute)
        return x.group(0).strip('strain: ').rstrip(' ').lower()
    else:
        print('Found {}'.format(sample_attribute))
        return np.nan


def summary_starlogs_over_runs(directory, list_of_srr):
    df = pd.DataFrame()
    files_not_found = []
    for run in list_of_srr:
        if not os.path.isfile(os.path.join(directory, run + '.tsv')):
            files_not_found.append(run)
            continue
        temp_df = pd.read_table(os.path.join(directory, run + '.tsv'))
        df = pd.concat([df, temp_df])
    return df, files_not_found


def get_enrichment_cds_stats(pickle_file):
    data = pickle.load(open(pickle_file, 'rb'))
    mean = np.nanmean(data.values())
    median = np.nanmedian(data.values())
    stddev = np.nanstd(data.values())
    minn = np.nanmin(data.values())
    maxx = np.nanmax(data.values())
    return minx, maxx, mean, median, stddev


def get_fragment_enrichment_score(txt_file):
    with open(txt_file) as fh:
        data = fh.read()
    enrichment = data.strip('\(').strip('\)').strip(' ').strip()
    enrichment, pval = enrichment.split(',')
    if 'nan' not in enrichment:
        return float(enrichment.strip('Enrichment: ')), float(
            pval.strip(')').strip('pval: '))
    else:
        return np.nan, 1
