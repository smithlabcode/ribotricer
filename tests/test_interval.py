"""Tests for the interval module."""

import pytest

from ribotricer.interval import Interval


class TestInterval:
    """Tests for the Interval class."""

    def test_interval_creation(self):
        """Test basic interval creation."""
        interval = Interval("chr1", 100, 200, "+")
        assert interval.chrom == "chr1"
        assert interval.start == 100
        assert interval.end == 200
        assert interval.strand == "+"

    def test_interval_length(self):
        """Test interval length calculation."""
        interval = Interval("chr1", 100, 200, "+")
        assert len(interval) == 101  # 1-based, closed interval

    def test_interval_equality(self):
        """Test interval equality comparison."""
        int1 = Interval("chr1", 100, 200, "+")
        int2 = Interval("chr1", 100, 200, "+")
        int3 = Interval("chr1", 100, 200, "-")

        assert int1 == int2
        assert int1 != int3

    def test_interval_str_representation(self):
        """Test string representation of interval."""
        interval = Interval("chr1", 100, 200, "+")
        str_repr = str(interval)
        assert "chr1" in str_repr
        assert "100" in str_repr
        assert "200" in str_repr

    def test_interval_hash(self):
        """Test that intervals can be hashed (for use in sets/dicts)."""
        int1 = Interval("chr1", 100, 200, "+")
        int2 = Interval("chr1", 100, 200, "+")

        # Same intervals should have same hash
        assert hash(int1) == hash(int2)

        # Should be usable in a set
        interval_set = {int1, int2}
        assert len(interval_set) == 1
