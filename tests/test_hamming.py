"""Tests for Hamming distance computation.

Tests trajectory.hamming_distance and branches.calculate_hamming_distance
handle edge cases correctly and agree with each other.
"""

from branches import calculate_hamming_distance
from trajectory import hamming_distance


class TestHammingDistance:
    """Test Hamming distance implementations and edge cases."""

    def test_identical_sequences(self):
        assert hamming_distance("ATCGATCGAT", "ATCGATCGAT") == 0

    def test_worked_example_x_y(self, sequences):
        """X vs Y: 1 difference at position 4 (G→A)."""
        assert hamming_distance(sequences["X"], sequences["Y"]) == 1

    def test_worked_example_y_a(self, sequences):
        """Y vs A: 2 differences (pos 4: A→G reversion, pos 6: T→G)."""
        assert hamming_distance(sequences["Y"], sequences["A"]) == 2

    def test_worked_example_y_b(self, sequences):
        """Y vs B: 1 difference (pos 9: A→G)."""
        assert hamming_distance(sequences["Y"], sequences["B"]) == 1

    def test_worked_example_x_c(self, sequences):
        """X vs C: 3 differences (pos 2: T→G, pos 5: A→G, pos 10: T→C)."""
        assert hamming_distance(sequences["X"], sequences["C"]) == 3

    def test_gaps_ignored(self):
        """Positions with gaps should not count as differences."""
        assert hamming_distance("ATCG-T", "ATCG-T") == 0
        assert hamming_distance("ATCG-T", "ATCGAT") == 0

    def test_ambiguous_n_ignored(self):
        """Positions with N should not count as differences."""
        assert hamming_distance("ATCGNT", "ATCGAT") == 0

    def test_mixed_gaps_and_diffs(self):
        """Only ACGT-vs-ACGT differences should be counted."""
        # One real diff (pos 1: A→G) + gap position + N position
        assert hamming_distance("A-NT", "G-AT") == 1

    def test_both_functions_agree(self, sequences):
        """trajectory.hamming_distance and branches.calculate_hamming_distance agree."""
        nodes = list(sequences.keys())
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                s1 = sequences[nodes[i]]
                s2 = sequences[nodes[j]]
                assert calculate_hamming_distance(s1, s2) == hamming_distance(s1, s2), \
                    f"Mismatch for {nodes[i]} vs {nodes[j]}"
