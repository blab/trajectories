"""Shared fixtures for trajectory pipeline tests.

Fixtures are built from the worked example in notes/data_format.md:
a tree with 3 tips (A, B, C), root X, internal node Y, and 10-nt sequences.
"""

import sys
import os

# Add scripts directory to path so we can import modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))

import pytest


@pytest.fixture
def sequences():
    """Sequences from the worked example in notes/data_format.md."""
    return {
        "X": "ATCGATCGAT",
        "Y": "ATCAATCGAT",
        "A": "ATCAAGCAAT",
        "B": "ATCAATCGGT",
        "C": "AGCGGTCGAC",
    }


@pytest.fixture
def parent_of():
    """Parent-child relationships from the worked example."""
    return {"Y": "X", "A": "Y", "B": "Y", "C": "X"}


@pytest.fixture
def hamming_of():
    """Branch Hamming distances from the worked example."""
    return {
        ("X", "Y"): 1,
        ("Y", "A"): 2,
        ("Y", "B"): 1,
        ("X", "C"): 3,
    }


@pytest.fixture
def gap_sequences():
    """Sequences with gaps to test gap-masking behavior.

    Z is a zero-distance node between Y and A with a different gap pattern.
    Y has a gap at position 9; Z fills it in. This means:
      _hamming_no_gaps(Y, A) = 2  (pos 9 skipped because gap in Y)
      _hamming_no_gaps(Z, A) = 3  (pos 9 counted because no gap in Z)
    The trajectory builder bug used Z→A (3) instead of Y→A (2) after skipping Z.
    """
    return {
        "X": "ATCGATCG-T",
        "Y": "ATCAATCG-T",  # 1 ACGT diff from X (pos 4: G→A)
        "Z": "ATCAATCGAT",  # 0 ACGT diff from Y (pos 9 is gap in Y → skipped)
        "A": "ATCAAGCAGT",  # 3 ACGT diffs from Z (pos 6,8,9), 2 from Y (pos 6,8)
    }


@pytest.fixture
def gap_parent_of():
    """Parent-child relationships for the gap test tree."""
    return {"Y": "X", "Z": "Y", "A": "Z"}


@pytest.fixture
def gap_hamming_of():
    """Branch Hamming distances for the gap test tree.

    These are the tree-edge values (computed between direct parent-child pairs).
    """
    return {
        ("X", "Y"): 1,
        ("Y", "Z"): 0,
        ("Z", "A"): 3,
    }
