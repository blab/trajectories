"""Tests for the forwards trajectory builder.

The most critical invariant: for every consecutive pair of emitted frames,
the header N must equal _hamming_no_gaps(prev_emitted_seq, curr_emitted_seq).
A test asserting this property would have caught the trajectory builder bug
where tree-edge Hamming values were used instead of direct sequence comparison.
"""

from trajectory import build_trajectory_content, _hamming_no_gaps


def parse_fasta_headers(content):
    """Parse FASTA content and return list of (name, distance) tuples."""
    headers = []
    for line in content.splitlines():
        if line.startswith(">"):
            name, dist = line[1:].split("|")
            headers.append((name, int(dist)))
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
            name, _ = line[1:].split("|")
            current_name = name
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
        actual_dist = _hamming_no_gaps(prev_seq, curr_seq)
        assert header_dist == actual_dist, (
            f"Header N invariant violated: {curr_name} header says {header_dist}, "
            f"but _hamming_no_gaps({prev_name}, {curr_name}) = {actual_dist}"
        )


class TestBasicTrajectory:
    """Test trajectory builder against the worked example from data_format.md."""

    def test_path_x_y_a(self, sequences, hamming_of):
        """Path X→Y→A produces headers X|0, Y|1, A|2."""
        path = ["X", "Y", "A"]
        content, tip_dist, depth = build_trajectory_content(path, sequences, hamming_of)
        headers = parse_fasta_headers(content)
        assert headers == [("X", 0), ("Y", 1), ("A", 2)]
        assert depth == 3
        assert tip_dist == 3  # cumulative: 0 + 1 + 2

    def test_path_x_y_b(self, sequences, hamming_of):
        """Path X→Y→B produces headers X|0, Y|1, B|1."""
        path = ["X", "Y", "B"]
        content, tip_dist, depth = build_trajectory_content(path, sequences, hamming_of)
        headers = parse_fasta_headers(content)
        assert headers == [("X", 0), ("Y", 1), ("B", 1)]
        assert depth == 3

    def test_path_x_c(self, sequences, hamming_of):
        """Path X→C produces headers X|0, C|3."""
        path = ["X", "C"]
        content, tip_dist, depth = build_trajectory_content(path, sequences, hamming_of)
        headers = parse_fasta_headers(content)
        assert headers == [("X", 0), ("C", 3)]
        assert depth == 2
        assert tip_dist == 3

    def test_header_invariant_all_paths(self, sequences, hamming_of):
        """Header N invariant holds for all paths in the worked example."""
        for path in [["X", "Y", "A"], ["X", "Y", "B"], ["X", "C"]]:
            content, _, _ = build_trajectory_content(path, sequences, hamming_of)
            assert_header_invariant(content)

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
        names = [name for name, _ in headers]
        assert "Z" not in names
        assert names == ["X", "Y", "A"]
        assert depth == 3

    def test_skip_preserves_header_invariant(self, sequences, hamming_of):
        """Header N invariant holds after skipping zero-distance nodes."""
        seqs = dict(sequences)
        seqs["Z"] = sequences["Y"]
        h = dict(hamming_of)
        h[("Y", "Z")] = 0
        h[("Z", "A")] = 2

        path = ["X", "Y", "Z", "A"]
        content, _, _ = build_trajectory_content(path, seqs, h)
        assert_header_invariant(content)


class TestZeroDistanceSkipWithGaps:
    """Regression test for the trajectory builder bug.

    When zero-distance node Z is skipped and Z has a different gap pattern
    than the last emitted node Y, the header for A must use
    _hamming_no_gaps(Y, A), not _hamming_no_gaps(Z, A). These differ because
    Z fills a gap that Y has at position 9, making position 9 visible in
    comparisons against Z but not against Y.
    """

    def test_gap_change_in_skipped_node(self, gap_sequences, gap_hamming_of):
        """Header N uses last emitted sequence, not the skipped node's sequence."""
        path = ["X", "Y", "Z", "A"]
        content, _, depth = build_trajectory_content(path, gap_sequences, gap_hamming_of)
        headers = parse_fasta_headers(content)

        # Z should be skipped
        names = [name for name, _ in headers]
        assert "Z" not in names
        assert names == ["X", "Y", "A"]

        # Header for A should be _hamming_no_gaps(Y, A) = 2, NOT tree-edge Z→A = 3
        a_header_dist = headers[2][1]
        expected = _hamming_no_gaps(gap_sequences["Y"], gap_sequences["A"])
        assert a_header_dist == expected == 2

    def test_header_invariant_with_gap_changes(self, gap_sequences, gap_hamming_of):
        """Header N invariant holds when skipped nodes have different gap patterns."""
        path = ["X", "Y", "Z", "A"]
        content, _, _ = build_trajectory_content(path, gap_sequences, gap_hamming_of)
        assert_header_invariant(content)


class TestCollapse:
    """Test collapse of zero-distance tip onto its parent."""

    def test_collapse_zero_distance_tip(self, sequences, hamming_of):
        """Tip T with zero distance from Y replaces Y in output."""
        seqs = dict(sequences)
        seqs["T"] = sequences["Y"]  # identical to Y
        h = dict(hamming_of)
        h[("Y", "T")] = 0

        path = ["X", "Y", "T"]
        content, _, depth = build_trajectory_content(path, seqs, h)
        headers = parse_fasta_headers(content)

        # Y should be replaced by T
        names = [name for name, _ in headers]
        assert "Y" not in names
        assert names == ["X", "T"]
        assert depth == 2

    def test_collapse_preserves_header_invariant(self, sequences, hamming_of):
        """Header N invariant holds after collapsing a tip."""
        seqs = dict(sequences)
        seqs["T"] = sequences["Y"]
        h = dict(hamming_of)
        h[("Y", "T")] = 0

        path = ["X", "Y", "T"]
        content, _, _ = build_trajectory_content(path, seqs, h)
        assert_header_invariant(content)

    def test_collapse_with_gap_changes(self):
        """Collapsed tip with different gap pattern gets correct header N."""
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

        # Y replaced by T; header for T = _hamming_no_gaps(X, T) = 1
        assert [name for name, _ in headers] == ["X", "T"]
        assert headers[1][1] == _hamming_no_gaps(seqs["X"], seqs["T"])
        assert_header_invariant(content)


class TestEdgeCases:
    """Test edge cases for trajectory builder."""

    def test_single_node_path(self, sequences, hamming_of):
        """Single-node path (just root) produces one frame."""
        path = ["X"]
        content, tip_dist, depth = build_trajectory_content(path, sequences, hamming_of)
        headers = parse_fasta_headers(content)
        assert headers == [("X", 0)]
        assert depth == 1
        assert tip_dist == 0

    def test_two_node_path(self, sequences, hamming_of):
        """Two-node path (root + child) produces two frames."""
        path = ["X", "Y"]
        content, tip_dist, depth = build_trajectory_content(path, sequences, hamming_of)
        headers = parse_fasta_headers(content)
        assert headers == [("X", 0), ("Y", 1)]
        assert depth == 2

    def test_missing_sequence_skipped(self, sequences, hamming_of):
        """Node with missing sequence is skipped gracefully."""
        seqs = {k: v for k, v in sequences.items() if k != "Y"}
        path = ["X", "Y", "A"]
        content, _, depth = build_trajectory_content(path, seqs, hamming_of)
        headers = parse_fasta_headers(content)
        names = [name for name, _ in headers]
        assert "Y" not in names
        assert "X" in names
        assert "A" in names
        assert_header_invariant(content)
