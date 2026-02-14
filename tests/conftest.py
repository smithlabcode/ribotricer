"""Pytest configuration and fixtures for ribotricer tests."""

from __future__ import annotations

import pytest

from ribotricer.interval import Interval


@pytest.fixture
def sample_interval() -> Interval:
    """Create a sample interval for testing."""
    return Interval("chr1", 100, 200, "+")


@pytest.fixture
def sample_intervals() -> list[Interval]:
    """Create a list of sample intervals for testing."""
    return [
        Interval("chr1", 100, 200, "+"),
        Interval("chr1", 300, 400, "+"),
        Interval("chr1", 500, 600, "+"),
    ]
