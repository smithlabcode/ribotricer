from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os

__GENOMES_DB__ = ['hg38', 'mm10']

__TAXON_ID_MAP__ = {
    '6239': 'C.elegans',
    '9606': 'H.sapiens',
    '10090': 'M.musculus',
    '9955': 'D.rerio',
    '559292': 'S.cerevisiase S288C',
    '7227': 'D.melanogaster',
    '10116': 'R.norvegicus'
}


def _get_bed(bedname, genome='hg38'):
    """Load bed from annotation.
    """
    annotation_dir = os.path.join(
        os.path.dirname(__file__), 'annotation', genome)
    annotation_file = os.path.join(annotation_dir, bedname + '.bed.gz')
    return annotation_file


def _get_sizes(genome):
    """Load genome sizes file"""
    annotation_dir = os.path.join(
        os.path.dirname(__file__), 'annotation', genome)
    return os.path.join(annotation_dir, '{}.chrom.sizes'.format(genome))
