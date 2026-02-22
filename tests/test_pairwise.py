"""Tests for pairwise trajectory generation."""

from pairwise_trajectory import calculate_hamming_distance, generate_pairs, build_pairwise_content


class TestPairwiseHamming:
    """Test pairwise Hamming distance computation."""

    def test_a_vs_b(self, sequences):
        """A vs B: 3 differences (pos 6, 8, 9)."""
        assert calculate_hamming_distance(sequences["A"], sequences["B"]) == 3

    def test_a_vs_c(self, sequences):
        """A vs C: 6 differences (pos 2, 4, 5, 6, 8, 10)."""
        assert calculate_hamming_distance(sequences["A"], sequences["C"]) == 6

    def test_gaps_ignored(self):
        """Gaps are not counted as differences."""
        assert calculate_hamming_distance("ATCG-T", "ATCGAT") == 0
        assert calculate_hamming_distance("A-CG", "G-CG") == 1

    def test_n_ignored(self):
        """Ambiguous bases N are not counted as differences."""
        assert calculate_hamming_distance("ATCGNT", "ATCGAT") == 0


class TestGeneratePairs:
    """Test pair generation logic."""

    def test_all_pairs_count(self):
        """All pairs of 4 items = C(4,2) = 6 pairs."""
        pairs = generate_pairs(["A", "B", "C", "D"])
        assert len(pairs) == 6

    def test_limited_pairs(self):
        """Limit returns the requested number of pairs."""
        pairs = generate_pairs(["A", "B", "C", "D"], limit=3, seed=42)
        assert len(pairs) == 3

    def test_pairs_are_unique(self):
        """All generated pairs are distinct."""
        pairs = generate_pairs(["A", "B", "C", "D", "E"], limit=5, seed=42)
        assert len(pairs) == len(set(pairs))

    def test_deterministic_with_seed(self):
        """Same seed produces the same pairs."""
        pairs1 = generate_pairs(list(range(20)), limit=5, seed=42)
        pairs2 = generate_pairs(list(range(20)), limit=5, seed=42)
        assert pairs1 == pairs2

    def test_two_items(self):
        """Two items produce exactly one pair."""
        pairs = generate_pairs(["A", "B"])
        assert len(pairs) == 1
        assert pairs[0] == ("A", "B")


class TestBuildPairwiseContent:
    """Test pairwise FASTA content building."""

    def test_basic_content(self, sequences):
        """First tip gets |0, second gets |{hamming}."""
        content, hamming = build_pairwise_content("A", "B", sequences)
        assert hamming == 3
        assert content.startswith(">A|0\n")
        assert ">B|3\n" in content

    def test_header_matches_direct_hamming(self, sequences):
        """Header distance equals direct Hamming between tip sequences."""
        for t1, t2 in [("A", "B"), ("A", "C"), ("B", "C")]:
            content, hamming = build_pairwise_content(t1, t2, sequences)
            expected = calculate_hamming_distance(sequences[t1], sequences[t2])
            assert hamming == expected

    def test_missing_sequence_returns_none(self, sequences):
        """Missing sequence returns (None, None)."""
        content, hamming = build_pairwise_content("A", "missing", sequences)
        assert content is None
        assert hamming is None

    def test_sequences_in_output(self, sequences):
        """Output contains the actual tip sequences."""
        content, _ = build_pairwise_content("A", "B", sequences)
        # 10-char sequences won't be line-wrapped at 60
        assert sequences["A"] in content
        assert sequences["B"] in content
