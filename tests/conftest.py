"""Pytest fixtures for simulate_data tests."""

import sys
from pathlib import Path

import pytest

# Ensure src/ is on the path when running tests
SRC_DIR = Path(__file__).parent.parent / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


@pytest.fixture
def mini_genome_path():
    """Path to the mini genome test FASTA."""
    return Path(__file__).parent / "data" / "mini_genome.fa"


@pytest.fixture
def mini_te_path():
    """Path to the mini TE consensus test FASTA."""
    return Path(__file__).parent / "data" / "mini_te.fa"


@pytest.fixture
def mini_te_bed_path():
    """Path to the mini TE positions test BED file."""
    return Path(__file__).parent / "data" / "mini_te.bed"
