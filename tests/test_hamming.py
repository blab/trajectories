"""Tests for Hamming distance computation.

Tests both implementations: branches.calculate_hamming_distance and
trajectory._hamming_no_gaps. They should produce identical results.
"""

from branches import calculate_hamming_distance
from trajectory import _hamming_no_gaps


class TestHammingDistance:
    """Test both Hamming implementations agree and handle edge cases."""

    def test_identical_sequences(self):
        assert calculate_hamming_distance("ATCGATCGAT", "ATCGATCGAT") == 0
        assert _hamming_no_gaps("ATCGATCGAT", "ATCGATCGAT") == 0

    def test_worked_example_x_y(self, sequences):
        """X vs Y: 1 difference at position 4 (G→A)."""
        assert calculate_hamming_distance(sequences["X"], sequences["Y"]) == 1
        assert _hamming_no_gaps(sequences["X"], sequences["Y"]) == 1

    def test_worked_example_y_a(self, sequences):
        """Y vs A: 2 differences (pos 4: A→G reversion, pos 6: T→G)."""
        assert calculate_hamming_distance(sequences["Y"], sequences["A"]) == 2
        assert _hamming_no_gaps(sequences["Y"], sequences["A"]) == 2

    def test_worked_example_y_b(self, sequences):
        """Y vs B: 1 difference (pos 9: A→G)."""
        assert calculate_hamming_distance(sequences["Y"], sequences["B"]) == 1
        assert _hamming_no_gaps(sequences["Y"], sequences["B"]) == 1

    def test_worked_example_x_c(self, sequences):
        """X vs C: 3 differences (pos 2: T→G, pos 5: A→G, pos 10: T→C)."""
        assert calculate_hamming_distance(sequences["X"], sequences["C"]) == 3
        assert _hamming_no_gaps(sequences["X"], sequences["C"]) == 3

    def test_gaps_ignored(self):
        """Positions with gaps should not count as differences."""
        assert calculate_hamming_distance("ATCG-T", "ATCG-T") == 0
        assert calculate_hamming_distance("ATCG-T", "ATCGAT") == 0
        assert _hamming_no_gaps("ATCG-T", "ATCGAT") == 0

    def test_ambiguous_n_ignored(self):
        """Positions with N should not count as differences."""
        assert calculate_hamming_distance("ATCGNT", "ATCGAT") == 0
        assert _hamming_no_gaps("ATCGNT", "ATCGAT") == 0

    def test_mixed_gaps_and_diffs(self):
        """Only ACGT-vs-ACGT differences should be counted."""
        # One real diff (pos 1: A→G) + gap position + N position
        assert calculate_hamming_distance("A-NT", "G-AT") == 1
        assert _hamming_no_gaps("A-NT", "G-AT") == 1

    def test_both_functions_agree(self, sequences):
        """Both Hamming implementations produce identical results on all pairs."""
        nodes = list(sequences.keys())
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                s1 = sequences[nodes[i]]
                s2 = sequences[nodes[j]]
                assert calculate_hamming_distance(s1, s2) == _hamming_no_gaps(s1, s2), \
                    f"Mismatch for {nodes[i]} vs {nodes[j]}"
