"""Tests for the ORF module."""

from __future__ import annotations

from ribotricer.interval import Interval
from ribotricer.orf import ORF


class TestORF:
    """Tests for the ORF class."""

    def test_orf_creation(self) -> None:
        """Test basic ORF creation."""
        intervals = [
            Interval("chr1", 100, 200, "+"),
            Interval("chr1", 300, 400, "+"),
        ]
        orf = ORF(
            category="annotated",
            transcript_id="tx1",
            transcript_type="protein_coding",
            gene_id="gene1",
            gene_name="Gene1",
            gene_type="protein_coding",
            chrom="chr1",
            strand="+",
            intervals=intervals,
            seq="ATGAAA",
        )

        assert orf.tid == "tx1"
        assert orf.gid == "gene1"
        assert orf.category == "annotated"
        assert orf.start_codon == "ATG"
        assert len(orf.intervals) == 2

    def test_orf_strand(self) -> None:
        """Test ORF strand property."""
        intervals = [Interval("chr1", 100, 200, "+")]
        orf = ORF(
            category="annotated",
            transcript_id="tx1",
            transcript_type="protein_coding",
            gene_id="gene1",
            gene_name="Gene1",
            gene_type="protein_coding",
            chrom="chr1",
            strand="+",
            intervals=intervals,
            seq="ATG",
        )

        assert orf.strand == "+"

    def test_orf_chrom(self) -> None:
        """Test ORF chromosome property."""
        intervals = [Interval("chr1", 100, 200, "+")]
        orf = ORF(
            category="annotated",
            transcript_id="tx1",
            transcript_type="protein_coding",
            gene_id="gene1",
            gene_name="Gene1",
            gene_type="protein_coding",
            chrom="chr1",
            strand="+",
            intervals=intervals,
            seq="ATG",
        )

        assert orf.chrom == "chr1"

    def test_orf_from_string(self) -> None:
        """Test ORF parsing from index file line."""
        line = (
            "tx1_100_200_101\tannotated\ttx1\tprotein_coding\t"
            "gene1\tGene1\tprotein_coding\tchr1\t+\tATG\t100-200"
        )
        orf = ORF.from_string(line)

        assert orf is not None
        assert orf.tid == "tx1"
        assert orf.category == "annotated"
        assert orf.start_codon == "ATG"
        assert orf.chrom == "chr1"
        assert orf.strand == "+"

    def test_orf_start_codon_short_seq(self) -> None:
        """Test that start_codon returns None for short sequences."""
        intervals = [Interval("chr1", 100, 102, "+")]
        orf = ORF(
            category="annotated",
            transcript_id="tx1",
            transcript_type="protein_coding",
            gene_id="gene1",
            gene_name="Gene1",
            gene_type="protein_coding",
            chrom="chr1",
            strand="+",
            intervals=intervals,
            seq="AT",  # Only 2 bases
        )

        assert orf.start_codon is None
