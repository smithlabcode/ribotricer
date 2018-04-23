'''Utility function for generating bed from GTFs'''
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from collections import defaultdict
import copy
import logging
import re
import warnings

import gffutils
import pybedtools
import pandas as pd


def features_to_df(feature_list):
    features = map(
        lambda x: (x.chrom, x.start, x.end, x.strand, x.feature_type),
        feature_list)
    return features


def create_gene_dict(db):
    '''
    Store each feature line db.all_features() as a dict of dicts

    Parameters
    ----------
    db : gffutils.FeatureDB
        Input database

    Returns
    -------
    gene_dict : dict
                Dictionary like {'gene_id': {'gene': gene_feature, 'transcript_id1': {'feature_type': [feature1, feature2,..]}}}
    '''
    gene_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for line_no, feature in enumerate(db.all_features()):
        gene_ids = feature.attributes['gene_id']
        feature_type = feature.featuretype
        if feature_type == 'gene':
            if len(gene_ids) != 1:
                logging.warning('Found multiple gene_ids on line {} in gtf'.
                                format(line_no))
                break
            else:
                gene_id = gene_ids[0]
                gene_dict[gene_id]['gene'] = feature
        else:
            transcript_ids = feature.attributes['transcript_id']

            for gene_id in gene_ids:
                for transcript_id in transcript_ids:
                    gene_dict[gene_id][transcript_id][feature_type].append(
                        feature)
    return gene_dict


def get_gene_list(gene_dict):
    return list(set(gene_dict.keys()))


def get_UTR_regions(gene_dict, gene_id, transcript, cds):
    '''Get UTR only regions

    Parameters
    ----------
    gene_dict : dict
                Output of create_gene_dict()
    gene_id : str
              Gene id
    transcript : str
                 Transcript id
    cds : str
          List of CDS

    Returns
    -------
    utr5_regions : list<feature>
                   List of gffutils.feature representing 5'UTR regions

    utr3_regions : list<feature>
                   List of gffutils.feature representing 5'UTR regions

    '''
    if len(cds) == 0:
        return [], []
    utr5_regions = []
    utr3_regions = []
    utrs = gene_dict[gene_id][transcript]['UTR']
    first_cds = cds[0]
    last_cds = cds[-1]
    for utr in utrs:
        # Push all cds at once
        # Sort later to remove duplicates
        strand = utr.strand
        if strand == '+':
            if utr.stop < first_cds.start:
                utr.feature_type = 'five_prime_UTR'
                utr5_regions.append(utr)
            elif utr.start > last_cds.stop:
                utr.feature_type = 'three_prime_UTR'
                utr3_regions.append(utr)
            else:
                raise RuntimeError('Error with cds')
        elif strand == '-':
            if utr.stop < first_cds.start:
                utr.feature_type = 'three_prime_UTR'
                utr3_regions.append(utr)
            elif utr.start > last_cds.stop:
                utr.feature_type = 'five_prime_UTR'
                utr5_regions.append(utr)
            else:
                raise RuntimeError('Error with cds')
    return utr5_regions, utr3_regions


def create_bed(regions, bedtype='0'):
    '''Create bed from list of regions.

    Parameters
    ----------
    regions : array like
              List of regions
    bedtype : 0 or 1
              0-Based or 1-based coordinate of the BED
    '''
    bedstr = ''
    for region in regions:
        assert len(region.attributes['gene_id']) == 1
        # GTF start is 1-based, so shift by one while writing
        # to 0-based BED format
        if bedtype == '0':
            start = region.start - 1
        else:
            start = region.start
        bedstr += '{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            region.chrom, start, region.stop,
            re.sub('\.\d+', '',
                   region.attributes['gene_id'][0]), '.', region.strand)
    return bedstr


def rename_regions(regions, gene_id):
    '''Name all regions one gene_id

    Parameters
    ----------
    regions : list<gffutils.Feature>
              Input
    gene_id : string
              Input

    Returns
    -------
    rename_regions : list<gffutils.Feature>
    '''
    regions = list(regions)
    if len(regions) == 0:
        return []
    for region in regions:
        region.attributes['gene_id'] = gene_id
    return regions


def merge_regions(db, regions):
    '''Merge overlapping regions respecting strand

    Parameters
    ----------
    db : gffutils.FeatureDB
         Input database

    regions : list<gffutils.Feature>
              List of regions to merge
    '''
    if len(regions) == 0:
        return []
    merged = db.merge(sorted(list(regions), key=lambda x: x.start))
    return merged


def merge_regions_nostrand(db, regions):
    '''Merge overlapping regions ignoring strand

    Parameters
    ----------
    db : gffutils.FeatureDB
         Input database

    regions : list<gffutils.Feature>
              List of regions to merge
    '''
    if len(regions) == 0:
        return []
    merged = db.merge(
        sorted(list(regions), key=lambda x: x.start), ignore_strand=True)
    return merged


def create_all_bed(gtf_file, prefix):
    '''Create beds from GTF

    Parameters
    ----------
    gtf_file : string
               Path to gtf file
    prefix : string
             Output prefix

    Returns
    -------
    beds : all beds
    '''
    db = gffutils.create_db(
        gtf_file,
        dbfn='{}.gffutils.db'.format(gtf_file),
        merge_strategy='merge',
        disable_infer_genes=True,
        disable_infer_transcripts=True)
    gene_dict = create_gene_dict(db)
    utr5_bed = ''
    utr3_bed = ''
    gene_bed = ''
    exon_bed = ''
    intron_bed = ''
    start_codon_bed = ''
    stop_codon_bed = ''
    cds_bed = ''

    gene_list = []
    all_utr5 = []
    all_utr3 = []

    for gene_id in get_gene_list(gene_dict):
        gene_list.append(gene_dict[gene_id]['gene'])

        utr5_regions, utr3_regions = [], []
        exon_regions, intron_regions = [], []
        star_codon_regions, stop_codon_regions = [], []
        cds_regions = []

        for feature in gene_dict[gene_id].keys():
            if feature == 'gene':
                continue
            cds = list(gene_dict[gene_id][feature]['CDS'])
            exons = list(gene_dict[gene_id][feature]['exon'])
            merged_exons = merge_regions(db, exons)
            introns = db.interfeatures(merged_exons)
            utr5_region, utr3_region = get_UTR_regions(gene_dict, gene_id,
                                                       feature, cds)
            utr5_regions += utr5_region
            all_utr5 += utr5_region
            all_utr3 += utr3_region
            utr3_regions += utr3_region
            exon_regions += exons
            intron_regions += introns
            cds_regions += cds

        merged_utr5 = merge_regions(db, utr5_regions)
        renamed_utr5 = rename_regions(merged_utr5, gene_id)

        merged_utr3 = merge_regions(db, utr3_regions)
        renamed_utr3 = rename_regions(merged_utr3, gene_id)

        merged_exons = merge_regions(db, exon_regions)
        renamed_exons = rename_regions(merged_exons, gene_id)

        merged_introns = merge_regions(db, intron_regions)
        renamed_introns = rename_regions(merged_introns, gene_id)

        merged_cds = merge_regions(db, cds_regions)
        renamed_cds = rename_regions(merged_cds, gene_id)

        utr3_bed += create_bed(renamed_utr3)
        utr5_bed += create_bed(renamed_utr5)
        exon_bed += create_bed(renamed_exons)
        intron_bed += create_bed(renamed_introns)
        cds_bed += create_bed(renamed_cds)
    all_utrs = all_utr5 + all_utr3
    all_utr_list = features_to_df(all_utrs)
    df = pd.DataFrame(
        all_utr_list, columns=[
            'chrom',
            'start',
            'end',
            'strand',
            'name',
        ])
    df = df.drop_duplicates()
    df = df.set_index(['chrom', 'start', 'end', 'strand'])
    with open('{}.modified_UTRS.gtf'.format(prefix), 'w') as fout:
        for d in db.directives:
            fout.write('## {0}\n'.format(d))
        for feature in db.all_features():
            if feature.featuretype == 'UTR':
                chrom, start, end, strand = feature.chrom, int(
                    feature.start), int(feature.end), feature.strand
                try:
                    feature.featuretype = str(
                        df.loc[(chrom, start, end, strand)].name[0])
                except KeyError:
                    print('something went wrong with feature: {}'.format(
                        feature))
            fout.write(str(feature) + '\n')

    gene_bed = create_bed(gene_list)
    gene_bedtool = pybedtools.BedTool(gene_bed, from_string=True)
    utr5_bedtool = pybedtools.BedTool(utr5_bed, from_string=True)
    utr3_bedtool = pybedtools.BedTool(utr3_bed, from_string=True)
    exon_bedtool = pybedtools.BedTool(exon_bed, from_string=True)
    intron_bedtool = pybedtools.BedTool(intron_bed, from_string=True)
    cds_bedtool = pybedtools.BedTool(cds_bed, from_string=True)

    gene_bedtool.remove_invalid().sort().saveas('{}.genes.bed'.format(prefix))
    utr5_bedtool.remove_invalid().sort().saveas('{}.UTR5.bed'.format(prefix))
    utr3_bedtool.remove_invalid().sort().saveas('{}.UTR3.bed'.format(prefix))
    exon_bedtool.remove_invalid().sort().saveas('{}.exon.bed'.format(prefix))
    intron_bedtool.remove_invalid().sort().saveas(
        '{}.intron.bed'.format(prefix))
    cds_bedtool.remove_invalid().sort().saveas('{}.cds.bed'.format(prefix))

    for gene_id in get_gene_list(gene_dict):
        start_codons = []
        stop_codons = []
        for start_codon in db.children(gene_id, featuretype='start_codon'):
            # 1 -based stop
            # 0-based start handled while converting to bed
            start_codon.stop = start_codon.start
            start_codons.append(start_codon)
        for stop_codon in db.children(gene_id, featuretype='stop_codon'):
            stop_codon.start = stop_codon.stop
            stop_codon.stop = stop_codon.stop + 1
            stop_codons.append(stop_codon)
        merged_start_codons = merge_regions(db, start_codons)
        renamed_start_codons = rename_regions(merged_start_codons, gene_id)
        merged_stop_codons = merge_regions(db, stop_codons)
        renamed_stop_codons = rename_regions(merged_stop_codons, gene_id)

        start_codon_bed += create_bed(renamed_start_codons)
        stop_codon_bed += create_bed(renamed_stop_codons)

    start_codon_bedtool = pybedtools.BedTool(start_codon_bed, from_string=True)
    stop_codon_bedtool = pybedtools.BedTool(stop_codon_bed, from_string=True)
    start_codon_bedtool.remove_invalid().sort().saveas(
        '{}.start_codon.bed'.format(prefix))
    stop_codon_bedtool.remove_invalid().sort().saveas(
        '{}.stop_codon.bed'.format(prefix))

    polyA_sites_bed = ''
    tss_sites_bed = ''
    for gene_id in get_gene_list(gene_dict):
        tss_sites = []
        polyA_sites = []
        for transcript in db.children(gene_id, featuretype='transcript'):
            start_t = copy.deepcopy(transcript)
            stop_t = copy.deepcopy(transcript)
            start_t.stop = start_t.start + 1
            stop_t.start = stop_t.stop

            if transcript.strand == '-':
                start_t, stop_t = stop_t, start_t
            polyA_sites.append(start_t)
            tss_sites.append(stop_t)
        merged_polyA_sites = merge_regions(db, polyA_sites)
        renamed_polyA_sites = rename_regions(merged_polyA_sites, gene_id)
        merged_tss_sites = merge_regions(db, tss_sites)
        renamed_tss_sites = rename_regions(merged_tss_sites, gene_id)
        polyA_sites_bed += create_bed(renamed_polyA_sites)
        tss_sites_bed += create_bed(renamed_tss_sites)

    polyA_sites_bedtool = pybedtools.BedTool(polyA_sites_bed, from_string=True)
    tss_sites_bedtool = pybedtools.BedTool(tss_sites_bed, from_string=True)
    polyA_sites_bedtool.remove_invalid().sort().saveas(
        '{}.polyA_sites.bed'.format(prefix))
    tss_sites_bedtool.remove_invalid().sort().saveas(
        '{}.tss_sites.bed'.format(prefix))
    rRNA_sites = []
    rRNA_bed = ''
    for gene_id in get_gene_list(gene_dict):
        for transcript in db.children(gene_id, featuretype='transcript'):
            if 'rRNA' in transcript.attributes['gene_type']:
                rRNA_sites.append(transcript)
        #renamed_rRNA_sites = rename_regions(rRNA_sites, gene_id)
        rRNA_bed += create_bed(rRNA_sites)

    rRNA_sites_bedtool = pybedtools.BedTool(rRNA_bed, from_string=True)
    rRNA_sites_bedtool.remove_invalid().sort().saveas(
        '{}.rRNA_sites.bed'.format(prefix))

    tRNA_sites = []
    tRNA_bed = ''
    for gene_id in get_gene_list(gene_dict):
        for transcript in db.children(gene_id, featuretype='transcript'):
            if 'tRNA' in transcript.attributes['gene_type'] or 'Mt_tRNA' in transcript.attributes['gene_type']:
                tRNA_sites.append(transcript)
        #merged_tRNA_sites = merge_regions_nostrand(db, tRNA_sites)
        #renamed_tRNA_sites = rename_regions(merged_tRNA_sites, gene_id)
        tRNA_bed += create_bed(tRNA_sites)

    tRNA_sites_bedtool = pybedtools.BedTool(tRNA_bed, from_string=True)
    tRNA_sites_bedtool.remove_invalid().sort().saveas(
        '{}.tRNA_sites.bed'.format(prefix))
