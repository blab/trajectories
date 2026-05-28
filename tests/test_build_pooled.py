"""Unit tests for scripts/build_pooled.py.

Uses a local-filesystem fake S3 client to exercise the full rollup
pipeline against small in-memory tar.zst fixtures without touching the
network.
"""

import io
import json
import os
import shutil
import tarfile
import tempfile

import pytest
import zstandard as zstd

import build_pooled


def make_tarball(path, members):
    """Write a tar.zst archive at path containing the given (name, content_bytes) members."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tar_buf = io.BytesIO()
    with tarfile.open(fileobj=tar_buf, mode='w') as tar:
        for name, content in members:
            info = tarfile.TarInfo(name=name)
            info.size = len(content)
            tar.addfile(info, io.BytesIO(content))
    compressed = zstd.ZstdCompressor().compress(tar_buf.getvalue())
    with open(path, 'wb') as f:
        f.write(compressed)


def read_tarball(path):
    """Read a tar.zst archive at path; return list of (name, content_bytes)."""
    return list(build_pooled.iter_tarball_members(path))


class FakePaginator:
    def __init__(self, root, bucket):
        self._root = root
        self._bucket = bucket

    def paginate(self, Bucket, Prefix):
        assert Bucket == self._bucket
        bucket_root = os.path.join(self._root, Bucket)
        contents = []
        prefix_dir = os.path.join(bucket_root, Prefix)
        # Walk the whole bucket and emit keys that start with Prefix
        for dirpath, _, filenames in os.walk(bucket_root):
            rel = os.path.relpath(dirpath, bucket_root)
            for fn in filenames:
                key = os.path.normpath(os.path.join(rel, fn)) if rel != '.' else fn
                key = key.replace(os.sep, '/')
                if key.startswith(Prefix):
                    contents.append({'Key': key})
        yield {'Contents': sorted(contents, key=lambda c: c['Key'])}


class FakeS3Client:
    """Minimal boto3-style client backed by a local directory tree.

    Layout: <root>/<bucket>/<key> on disk.
    """

    def __init__(self, root):
        self._root = root

    def get_paginator(self, op):
        assert op == 'list_objects_v2'
        return _PaginatorFactory(self._root)

    def download_file(self, Bucket, Key, Filename):
        src = os.path.join(self._root, Bucket, Key)
        os.makedirs(os.path.dirname(Filename), exist_ok=True)
        shutil.copy(src, Filename)


class _PaginatorFactory:
    def __init__(self, root):
        self._root = root

    def paginate(self, Bucket, Prefix):
        return FakePaginator(self._root, Bucket).paginate(Bucket=Bucket, Prefix=Prefix)


@pytest.fixture
def fixture_bucket(tmp_path):
    """Create a fake S3 bucket on disk with a known layout and return (root, bucket, source_prefix)."""
    bucket = "test-bucket"
    source_prefix = "trellis-trajectories"
    root = tmp_path / "s3"
    base = root / bucket / source_prefix

    # Source A — non-holdout-style: forwards train + forwards test
    make_tarball(
        str(base / "trellis-18aa-AAEC" / "forwards-train-000.tar.zst"),
        members=[
            ("TIP_0001.fasta", b">NODE_0000|0|0\nATGCAT\n>TIP_0001|2|2\nATGAAT\n"),
            ("TIP_0002.fasta", b">NODE_0000|0|0\nATGCAT\n>TIP_0002|1|1\nATGCAA\n"),
        ],
    )
    make_tarball(
        str(base / "trellis-18aa-AAEC" / "forwards-test-000.tar.zst"),
        members=[
            ("TIP_0003.fasta", b">NODE_0000|0|0\nATGCAT\n>TIP_0003|3|3\nAAACAT\n"),
        ],
    )

    # Source B — non-holdout-style
    make_tarball(
        str(base / "trellis-18aa-AAFV" / "forwards-train-000.tar.zst"),
        members=[
            ("TIP_0001.fasta", b">NODE_0000|0|0\nATGCAT\n>TIP_0001|1|1\nATGCAA\n"),
        ],
    )
    make_tarball(
        str(base / "trellis-18aa-AAFV" / "forwards-test-000.tar.zst"),
        members=[
            ("TIP_0002.fasta", b">NODE_0000|0|0\nATGCAT\n>TIP_0002|2|2\nATAAAT\n"),
        ],
    )

    # Source C — phylo-style nested path with hyphenated subdir
    make_tarball(
        str(base / "trellis-18aa-AAEC-phylo" / "forwards-train-000.tar.zst"),
        members=[
            ("node_x.fasta", b">root|0|0\nATGCAT\n>node_x|1|1\nATGCAA\n"),
        ],
    )

    # Pairwise on one source
    make_tarball(
        str(base / "trellis-18aa-AAEC" / "pairwise-train-000.tar.zst"),
        members=[
            ("a__b.fasta", b">a|0|0\nATGCAT\n>b|5|5\nGGCCTT\n"),
        ],
    )

    # Spurious existing pooled — should be EXCLUDED
    make_tarball(
        str(base / "pooled" / "forwards-train-000.tar.zst"),
        members=[
            ("ghost.fasta", b">ghost|0|0\nXXXX\n"),
        ],
    )

    # Non-shard junk that should be ignored (skipped, not errored)
    make_tarball(
        str(base / "trellis-18aa-AAEC" / "extra-train-000.tar.zst"),
        members=[("junk.fasta", b">junk|0|0\nNNNN\n")],
    )

    # Non-tarball file that should be ignored entirely
    (base / "trellis-18aa-AAEC" / "README.txt").write_text("ignore me")

    return str(root), bucket, source_prefix


def test_classify_key():
    assert build_pooled.classify_key(
        "trellis-trajectories/trellis-18aa-AAEC/forwards-train-000.tar.zst",
        "trellis-trajectories",
    ) == ("forwards", "train", "trellis-18aa-AAEC", "forwards-train-000.tar.zst")

    # Nested subpath should be flattened with hyphens
    assert build_pooled.classify_key(
        "trajectories/rdrp-flaviviridae-xs/subtrees/clade-A/forwards-test-000.tar.zst",
        "trajectories",
    ) == ("forwards", "test", "rdrp-flaviviridae-xs-subtrees-clade-A",
          "forwards-test-000.tar.zst")

    # Non-shard filename should return None
    assert build_pooled.classify_key(
        "trellis-trajectories/trellis-18aa-AAEC/extra-train-000.tar.zst",
        "trellis-trajectories",
    ) is None

    # Top-level key (no subpath) should return None
    assert build_pooled.classify_key(
        "trellis-trajectories/forwards-train-000.tar.zst",
        "trellis-trajectories",
    ) is None

    # Key outside source prefix should return None
    assert build_pooled.classify_key(
        "other-prefix/trellis-18aa-AAEC/forwards-train-000.tar.zst",
        "trellis-trajectories",
    ) is None


def test_count_fasta_frames_bases():
    content = ">h1|0|0\nATGCAT\n>h2|3|3\nAAAA\nCC\n"
    frames, bases = build_pooled.count_fasta_frames_bases(content)
    assert frames == 2
    assert bases == 6 + 6  # ATGCAT + AAAA + CC


def test_iter_tarball_members(tmp_path):
    path = tmp_path / "test.tar.zst"
    members = [("a.fasta", b"hello"), ("b.fasta", b"world")]
    make_tarball(str(path), members)
    result = list(build_pooled.iter_tarball_members(str(path)))
    assert result == members


def test_build_pooled_end_to_end(fixture_bucket, tmp_path):
    root, bucket, source_prefix = fixture_bucket
    s3 = FakeS3Client(root)
    output_dir = tmp_path / "pooled"

    metadata = build_pooled.build_pooled(
        s3=s3,
        bucket=bucket,
        source_prefix=source_prefix,
        output_dir=str(output_dir),
        shard_size=1000,
        staging_dir=str(tmp_path / "staging"),
        keep_staging=False,
        shuffle=False,
    )

    # forwards-train: AAEC(2) + AAFV(1) + AAEC-phylo(1) = 4 trajectories
    # forwards-test: AAEC(1) + AAFV(1) = 2 trajectories
    # pairwise-train: AAEC(1) = 1 trajectory
    splits = metadata['splits']
    assert splits['forwards-train']['trajectories'] == 4
    assert splits['forwards-test']['trajectories'] == 2
    assert splits['pairwise-train']['trajectories'] == 1
    assert 'pairwise-test' not in splits

    # Shard counts (small fixture, all fit in shard 0)
    assert splits['forwards-train']['shards'] == 1
    assert splits['forwards-test']['shards'] == 1
    assert splits['pairwise-train']['shards'] == 1

    # source_keys excludes the pooled/ ghost and the non-shard extra-train file
    keys = metadata['source_keys']
    assert all('pooled/' not in k for k in keys)
    assert all('extra-train' not in k for k in keys)
    # Should include the 6 legitimate shards
    assert len(keys) == 6

    # Pooled tar member names are prefixed and content is byte-for-byte
    train_shard = output_dir / "forwards-train-000.tar.zst"
    train_members = dict(read_tarball(str(train_shard)))
    assert set(train_members.keys()) == {
        "trellis-18aa-AAEC-TIP_0001.fasta",
        "trellis-18aa-AAEC-TIP_0002.fasta",
        "trellis-18aa-AAFV-TIP_0001.fasta",
        "trellis-18aa-AAEC-phylo-node_x.fasta",
    }
    # Spot-check byte-for-byte content preservation for one file
    assert train_members["trellis-18aa-AAEC-TIP_0001.fasta"] == (
        b">NODE_0000|0|0\nATGCAT\n>TIP_0001|2|2\nATGAAT\n"
    )

    # Phylo content preserved
    assert train_members["trellis-18aa-AAEC-phylo-node_x.fasta"] == (
        b">root|0|0\nATGCAT\n>node_x|1|1\nATGCAA\n"
    )

    # metadata.json was written
    with open(output_dir / "metadata.json") as f:
        loaded = json.load(f)
    assert loaded == metadata


def test_existing_pooled_is_excluded(fixture_bucket, tmp_path):
    root, bucket, source_prefix = fixture_bucket
    s3 = FakeS3Client(root)
    keys = build_pooled.list_source_keys(s3, bucket, source_prefix)
    assert all('pooled/' not in k for k in keys), keys


def test_byte_for_byte_fasta_preservation(fixture_bucket, tmp_path):
    """Critical: pooled tarballs must preserve FASTA header + sequence bytes exactly.
    No ligand info should leak into headers."""
    root, bucket, source_prefix = fixture_bucket
    s3 = FakeS3Client(root)
    output_dir = tmp_path / "pooled"

    build_pooled.build_pooled(
        s3=s3,
        bucket=bucket,
        source_prefix=source_prefix,
        output_dir=str(output_dir),
        shard_size=1000,
        staging_dir=str(tmp_path / "staging"),
        keep_staging=False,
        shuffle=False,
    )

    for shard_name in ("forwards-train-000.tar.zst", "forwards-test-000.tar.zst",
                       "pairwise-train-000.tar.zst"):
        members = read_tarball(str(output_dir / shard_name))
        for name, content in members:
            # Source identifier (ligand) must appear in tar member name (provenance)
            assert any(prefix in name for prefix in (
                "trellis-18aa-AAEC", "trellis-18aa-AAFV", "trellis-18aa-AAEC-phylo")), name
            # But must NOT appear inside the FASTA file (no training contamination)
            text = content.decode('utf-8')
            for header_line in (line for line in text.splitlines() if line.startswith('>')):
                assert "trellis" not in header_line, (
                    f"FASTA header leaked ligand info: {header_line!r}"
                )
                assert "AAEC" not in header_line, (
                    f"FASTA header leaked ligand info: {header_line!r}"
                )
                assert "AAFV" not in header_line, (
                    f"FASTA header leaked ligand info: {header_line!r}"
                )
