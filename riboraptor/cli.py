from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import sys
import click
click.disable_unicode_literals_warning = True
import six

from riboraptor.utils import bedgraph_to_bigwig
from riboraptor.utils import fragment_enrichment
from riboraptor.utils import gene_coverage
from riboraptor.utils import htseq_to_cpm
from riboraptor.utils import mapping_reads_summary
from riboraptor.utils import read_length_distribution

from riboraptor.plotting import plot_read_counts
from riboraptor.plotting import plot_fragment_dist
from click_help_colors import HelpColorsGroup, HelpColorsCommand

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(cls=HelpColorsGroup,
             help_headers_color='yellow',
             help_options_color='green')
@click.version_option(version='0.1.0')
def cli():
    """riboraptor: Tool for ribosome profiling analysis"""
    pass


@cli.command('rld', context_settings=CONTEXT_SETTINGS)
@click.option('--bam',
              help='Path to BAM file',
              required=True)
def rld_cmd(bam):
    counts = read_length_distribution(bam)
    for l, count in six.iteritems(dict(counts)):
        print('{}\t{}'.format(l, count))


@cli.command('fragment-enrichment', context_settings=CONTEXT_SETTINGS)
@click.option('--lrange',
              help='Fragments lengths to use for enrichment',
              default='28-32')
def fragment_enrichment_cmd(lrange):
    enrichment, pvalue = fragment_enrichment(sys.stdin.readlines(),
                                             lrange,
                                             True)
    print('({}, {})'.format(enrichment, pvalue))
