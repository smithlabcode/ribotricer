from __future__ import absolute_import
from __future__ import print_function
import pytest
import pybedtools

from riboraptor.genome import _get_bed

from riboraptor.count import gene_coverage
from riboraptor.count import count_uniq_mapping_reads
from riboraptor.fasta import FastaReader
from riboraptor.interval import Interval
from riboraptor.sequence import gene_sequence
from riboraptor.sequence import export_gene_sequences
from riboraptor.wig import WigReader


def test_interval():
    iv = Interval()
    assert (iv.chrom == None and iv.start == 0 and iv.end == 0
            and iv.strand == None)
    iv = Interval('chr1', 1, 15, '+')
    assert (iv.chrom == 'chr1' and iv.start == 1 and iv.end == 15
            and iv.strand == '+')


def test_wig():
    bw = WigReader('tests/data/SRX2536403_subsampled.unique.bigWig')
    intervals = [
        Interval('chr1', x[0], x[1], '-')
        for x in [(999594, 999595), (999597, 999598)]
    ]
    coverages = bw.query(intervals)
    assert (coverages[0] == [1] and coverages[1] == [1])


def test_fasta():
    fasta = FastaReader('tests/data/hg38.fa')
    chrom_length = fasta.chromosomes['chr1']
    assert (chrom_length == 89950)
    sequences = fasta.query(
        [Interval('chr1', 89948, 89950),
         Interval('chr1', 99, 100)])
    assert (sequences == ['TTA', 'AC'])

    assert (fasta.complement('TCGA') == 'AGCT')
    assert (fasta.reverse_complement('TCGA') == 'TCGA')


def test_gene_coverage():
    gene_name = 'ENSG00000035115'
    bed = 'riboraptor/annotation/hg38/cds.bed.gz'
    bw = 'tests/data/SRX2536403_subsampled.unique.bigWig'
    bed_df = pybedtools.BedTool(bed).sort().to_dataframe()
    bed_df['chrom'] = bed_df['chrom'].astype(str)
    bed_df['name'] = bed_df['name'].astype(str)
    bed_grouped = bed_df.groupby('name')
    gene_group = bed_df[bed_df['name'] == gene_name]

    coverage, offset_5p, offset_3p = gene_coverage(gene_group, bw, 30, 30)
    assert (coverage.sum() == 1)


def test_gene_coverage_from_internal_bed():
    gene_name = 'ENSG00000035115'
    bed = 'hg38_cds'
    bw = 'tests/data/SRX2536403_subsampled.unique.bigWig'

    genome, region_type = bed.lower().split('_')
    bed = _get_bed(region_type, genome)
    bed_df = pybedtools.BedTool(bed).sort().to_dataframe()
    bed_df['chrom'] = bed_df['chrom'].astype(str)
    bed_df['name'] = bed_df['name'].astype(str)
    bed_grouped = bed_df.groupby('name')
    gene_group = bed_df[bed_df['name'] == gene_name]

    coverage, offset_5p, offset_3p = gene_coverage(gene_group, bw, 30, 30)
    assert (coverage.sum() == 1)


def test_uniq_mapping_count():
    bam = 'tests/data/SRX2536403_subsampled.unique.bam'
    counts = count_uniq_mapping_reads(bam)
    assert (int(counts) == 358514)


if __name__ == '__main__':
    pytest.main([__file__])
