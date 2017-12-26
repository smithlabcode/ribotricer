import os

__GENOMES_DB__ = ['hg38', 'mm10']


def _load_bed(bedname, genome='hg38'):
    """Load bed from annotation.
    """
    annotation_dir = os.path.join(os.path.dirname(__file__),
                                  'annotation',
                                  genome)
    annotation_file = os.path.join(annotation_dir,
                                   bedname + '.bed.gz')
    return annotation_file


def _get_sizes(genome):
    """Load genome sizes file"""
    annotation_dir = os.path.join(os.path.dirname(__file__),
                                  'annotation',
                                  genome)
    return os.path.join(annotation_dir,
                        '{}.chrom.sizes'.format(genome))
