"""Tests for the forwards trajectory builder.

The most critical invariant: for every consecutive pair of emitted frames,
the header N must equal hamming_distance(prev_emitted_seq, curr_emitted_seq).
A test asserting this property would have caught the trajectory builder bug
where tree-edge Hamming values were used instead of direct sequence comparison.
"""

from trajectory import build_trajectory_content, hamming_distance


def parse_fasta_headers(content):
    """Parse FASTA content and return list of (name, branch_dist, direct_dist) tuples."""
    headers = []
    for line in content.splitlines():
        if line.startswith(">"):
            parts = line[1:].split("|")
            headers.append((parts[0], int(parts[1]), int(parts[2])))
    return headers


def parse_fasta_sequences(content):
    """Parse FASTA content and return list of (name, sequence) tuples."""
    entries = []
    current_name = None
    current_seq = []
    for line in content.splitlines():
        if line.startswith(">"):
            if current_name is not None:
                entries.append((current_name, "".join(current_seq)))
            parts = line[1:].split("|")
            current_name = parts[0]
            current_seq = []
        else:
            current_seq.append(line)
    if current_name is not None:
        entries.append((current_name, "".join(current_seq)))
    return entries


def assert_header_invariant(content):
    """Assert header N equals ACGT-only Hamming between consecutive emitted sequences.

    This is the single most important invariant in the pipeline.
    """
    entries = parse_fasta_sequences(content)
    headers = parse_fasta_headers(content)
    for i in range(1, len(entries)):
        prev_name, prev_seq = entries[i - 1]
        curr_name, curr_seq = entries[i]
        header_dist = headers[i][1]
        actual_dist = hamming_distance(prev_seq, curr_seq)
        assert header_dist == actual_dist, (
            f"Header N invariant violated: {curr_name} header says {header_dist}, "
            f"but hamming_distance({prev_name}, {curr_name}) = {actual_dist}"
        )


def assert_direct_distance_invariant(content):
    """Assert direct_dist equals hamming_distance(start_seq, current_seq) for every frame."""
    entries = parse_fasta_sequences(content)
    headers = parse_fasta_headers(content)
    if not entries:
        return
    start_name, start_seq = entries[0]
    for i, (curr_name, curr_seq) in enumerate(entries):
        direct_dist = headers[i][2]
        actual_dist = hamming_distance(start_seq, curr_seq)
        assert direct_dist == actual_dist, (
            f"Direct distance invariant violated: {curr_name} header says {direct_dist}, "
            f"but hamming_distance({start_name}, {curr_name}) = {actual_dist}"
        )


class TestBasicTrajectory:
    """Test trajectory builder against the worked example from data_format.md."""

    def test_path_x_y_a(self, sequences, hamming_of):
        """Path X→Y→A produces headers X|0|0, Y|1|1, A|2|1.

        A reverts position 4 (A→G, back to root X) and mutates position 6 (T→G),
        so branch distance from Y is 2 but direct distance from X is only 1.
        """
        path = ["X", "Y", "A"]
        content, tip_dist, depth = build_trajectory_content(path, sequences, hamming_of)
        headers = parse_fasta_headers(content)
        assert headers == [("X", 0, 0), ("Y", 1, 1), ("A", 2, 1)]
        assert depth == 3
        assert tip_dist == 3  # cumulative: 0 + 1 + 2

    def test_path_x_y_b(self, sequences, hamming_of):
        """Path X→Y→B produces headers X|0|0, Y|1|1, B|1|2."""
        path = ["X", "Y", "B"]
        content, tip_dist, depth = build_trajectory_content(path, sequences, hamming_of)
        headers = parse_fasta_headers(content)
        assert headers == [("X", 0, 0), ("Y", 1, 1), ("B", 1, 2)]
        assert depth == 3

    def test_path_x_c(self, sequences, hamming_of):
        """Path X→C produces headers X|0|0, C|3|3."""
        path = ["X", "C"]
        content, tip_dist, depth = build_trajectory_content(path, sequences, hamming_of)
        headers = parse_fasta_headers(content)
        assert headers == [("X", 0, 0), ("C", 3, 3)]
        assert depth == 2
        assert tip_dist == 3

    def test_header_invariant_all_paths(self, sequences, hamming_of):
        """Header invariants hold for all paths in the worked example."""
        for path in [["X", "Y", "A"], ["X", "Y", "B"], ["X", "C"]]:
            content, _, _ = build_trajectory_content(path, sequences, hamming_of)
            assert_header_invariant(content)
            assert_direct_distance_invariant(content)

    def test_sequences_in_output(self, sequences, hamming_of):
        """Emitted sequences match the input sequences."""
        path = ["X", "Y", "A"]
        content, _, _ = build_trajectory_content(path, sequences, hamming_of)
        entries = parse_fasta_sequences(content)
        for name, seq in entries:
            assert seq == sequences[name]


class TestZeroDistanceSkip:
    """Test skipping of zero-distance intermediate nodes."""

    def test_skip_zero_distance_node(self, sequences, hamming_of):
        """Zero-distance node Z between Y and A is skipped in output."""
        seqs = dict(sequences)
        seqs["Z"] = sequences["Y"]  # identical sequence
        h = dict(hamming_of)
        h[("Y", "Z")] = 0
        h[("Z", "A")] = 2  # same as Y→A

        path = ["X", "Y", "Z", "A"]
        content, _, depth = build_trajectory_content(path, seqs, h)
        headers = parse_fasta_headers(content)
        names = [name for name, _, _ in headers]
        assert "Z" not in names
        assert names == ["X", "Y", "A"]
        assert depth == 3

    def test_skip_preserves_header_invariant(self, sequences, hamming_of):
        """Header invariants hold after skipping zero-distance nodes."""
        seqs = dict(sequences)
        seqs["Z"] = sequences["Y"]
        h = dict(hamming_of)
        h[("Y", "Z")] = 0
        h[("Z", "A")] = 2

        path = ["X", "Y", "Z", "A"]
        content, _, _ = build_trajectory_content(path, seqs, h)
        assert_header_invariant(content)
        assert_direct_distance_invariant(content)


class TestZeroDistanceSkipWithGaps:
    """Regression test for the trajectory builder bug.

    When zero-distance node Z is skipped and Z has a different gap pattern
    than the last emitted node Y, the header for A must use
    hamming_distance(Y, A), not hamming_distance(Z, A). These differ because
    Z fills a gap that Y has at position 9, making position 9 visible in
    comparisons against Z but not against Y.
    """

    def test_gap_change_in_skipped_node(self, gap_sequences, gap_hamming_of):
        """Header N uses last emitted sequence, not the skipped node's sequence."""
        path = ["X", "Y", "Z", "A"]
        content, _, depth = build_trajectory_content(path, gap_sequences, gap_hamming_of)
        headers = parse_fasta_headers(content)

        # Z should be skipped
        names = [name for name, _, _ in headers]
        assert "Z" not in names
        assert names == ["X", "Y", "A"]

        # Header for A should be hamming_distance(Y, A) = 2, NOT tree-edge Z→A = 3
        a_header_dist = headers[2][1]
        expected = hamming_distance(gap_sequences["Y"], gap_sequences["A"])
        assert a_header_dist == expected == 2

    def test_header_invariant_with_gap_changes(self, gap_sequences, gap_hamming_of):
        """Header invariants hold when skipped nodes have different gap patterns."""
        path = ["X", "Y", "Z", "A"]
        content, _, _ = build_trajectory_content(path, gap_sequences, gap_hamming_of)
        assert_header_invariant(content)
        assert_direct_distance_invariant(content)


class TestZeroDistanceTip:
    """Test that zero-distance tips are relabeled (not emitted as extra frame)."""

    def test_zero_distance_tip_relabeled(self, sequences, hamming_of):
        """Tip T with zero distance from Y relabels Y's frame."""
        seqs = dict(sequences)
        seqs["T"] = sequences["Y"]  # identical to Y
        h = dict(hamming_of)
        h[("Y", "T")] = 0

        path = ["X", "Y", "T"]
        content, _, depth = build_trajectory_content(path, seqs, h)
        headers = parse_fasta_headers(content)

        names = [name for name, _, _ in headers]
        assert names == ["X", "T"]
        assert depth == 2
        assert headers == [("X", 0, 0), ("T", 1, 1)]

    def test_zero_distance_tip_preserves_header_invariant(self, sequences, hamming_of):
        """Header invariants hold with zero-distance tip."""
        seqs = dict(sequences)
        seqs["T"] = sequences["Y"]
        h = dict(hamming_of)
        h[("Y", "T")] = 0

        path = ["X", "Y", "T"]
        content, _, _ = build_trajectory_content(path, seqs, h)
        assert_header_invariant(content)
        assert_direct_distance_invariant(content)

    def test_zero_distance_tip_with_gap_changes(self):
        """Zero-distance tip with different gap pattern relabels predecessor."""
        seqs = {
            "X": "ATCGATCGAT",
            "Y": "ATCAATCG-T",  # 1 ACGT diff from X, gap at pos 9
            "T": "ATCAATCGAT",  # 0 ACGT diff from Y (gap filled), but is tip
        }
        h = {
            ("X", "Y"): 1,
            ("Y", "T"): 0,
        }

        path = ["X", "Y", "T"]
        content, _, depth = build_trajectory_content(path, seqs, h)
        headers = parse_fasta_headers(content)

        names = [name for name, _, _ in headers]
        assert names == ["X", "T"]
        assert depth == 2
        # T relabels Y's frame, so header has Y's distances
        assert headers[1] == ("T", 1, 1)
        assert_header_invariant(content)
        assert_direct_distance_invariant(content)

    def test_chain_of_zero_distance_to_tip(self):
        """Chain of zero-distance nodes to tip: N1 relabeled as T.

        N2, N3, and T are all zero-distance from N1, so they are skipped.
        N1 is relabeled with the tip name T.
        """
        seqs = {
            "R": "ATCGATCGAT",
            "N1": "ATCAATCGAT",  # 1 diff from R
            "N2": "ATCAATCGAT",  # 0 diff from N1 (intermediate, skipped)
            "N3": "ATCAATCGAT",  # 0 diff from N2 (intermediate, skipped)
            "T":  "ATCAATCGAT",  # 0 diff from N3 (tip, skipped & relabeled)
        }
        h = {
            ("R", "N1"): 1,
            ("N1", "N2"): 0,
            ("N2", "N3"): 0,
            ("N3", "T"): 0,
        }

        path = ["R", "N1", "N2", "N3", "T"]
        content, _, depth = build_trajectory_content(path, seqs, h)
        headers = parse_fasta_headers(content)

        names = [name for name, _, _ in headers]
        assert names == ["R", "T"]
        assert depth == 2
        assert headers == [("R", 0, 0), ("T", 1, 1)]
        assert_header_invariant(content)
        assert_direct_distance_invariant(content)


class TestEdgeCases:
    """Test edge cases for trajectory builder."""

    def test_single_node_path(self, sequences, hamming_of):
        """Single-node path (just root) produces one frame."""
        path = ["X"]
        content, tip_dist, depth = build_trajectory_content(path, sequences, hamming_of)
        headers = parse_fasta_headers(content)
        assert headers == [("X", 0, 0)]
        assert depth == 1
        assert tip_dist == 0

    def test_two_node_path(self, sequences, hamming_of):
        """Two-node path (root + child) produces two frames."""
        path = ["X", "Y"]
        content, tip_dist, depth = build_trajectory_content(path, sequences, hamming_of)
        headers = parse_fasta_headers(content)
        assert headers == [("X", 0, 0), ("Y", 1, 1)]
        assert depth == 2

    def test_missing_sequence_skipped(self, sequences, hamming_of):
        """Node with missing sequence is skipped gracefully."""
        seqs = {k: v for k, v in sequences.items() if k != "Y"}
        path = ["X", "Y", "A"]
        content, _, depth = build_trajectory_content(path, seqs, hamming_of)
        headers = parse_fasta_headers(content)
        names = [name for name, _, _ in headers]
        assert "Y" not in names
        assert "X" in names
        assert "A" in names
        assert_header_invariant(content)
        assert_direct_distance_invariant(content)
