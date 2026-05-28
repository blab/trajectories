"""Build pooled trajectory shards from per-source granular tar.zst files on S3.

Walks an S3 prefix for granular trajectory shards (e.g.
`s3://bucket/prefix/source-X/forwards-train-NNN.tar.zst`), downloads them
locally, and repackages their contents into fixed-size pooled shards under
a single output directory. Train/test labels are taken from each source
shard's filename; no re-splitting is performed.

Pooled tar member names are prefixed with the source subpath (flattened
with hyphens) for provenance, e.g. `trellis-18aa-AAEC-TIP_0042.fasta`.
FASTA file contents are passed through byte-for-byte; headers are never
modified.
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
from collections import defaultdict
from datetime import datetime, timezone

import boto3
import zstandard as zstd
from tqdm import tqdm

from shard_writer import ShardWriter


SHARD_RE = re.compile(r'^(forwards|pairwise)-(train|test)-\d+\.tar\.zst$')
POOLED_DIR = 'pooled'


def get_git_commit():
    """Return short git commit hash, with '-dirty' suffix if working tree has changes."""
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL,
        ).decode().strip()
        try:
            subprocess.check_call(
                ["git", "diff", "--quiet", "HEAD"],
                stderr=subprocess.DEVNULL,
            )
        except subprocess.CalledProcessError:
            commit += "-dirty"
        return commit
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def list_source_keys(s3, bucket, source_prefix):
    """List all granular tar.zst keys under source_prefix, excluding pooled/.

    Returns a list of S3 keys (full key strings).
    """
    normalized_prefix = source_prefix.rstrip('/') + '/'
    pooled_prefix = normalized_prefix + POOLED_DIR + '/'
    paginator = s3.get_paginator('list_objects_v2')
    keys = []
    for page in paginator.paginate(Bucket=bucket, Prefix=normalized_prefix):
        for obj in page.get('Contents', []):
            key = obj['Key']
            if not key.endswith('.tar.zst'):
                continue
            if key.startswith(pooled_prefix):
                continue
            keys.append(key)
    return keys


def classify_key(key, source_prefix):
    """Parse a source key into (traj_type, split, flat_subpath, shard_filename).

    Returns None for keys whose filename doesn't match the expected pattern.
    """
    normalized_prefix = source_prefix.rstrip('/') + '/'
    if not key.startswith(normalized_prefix):
        return None
    relative = key[len(normalized_prefix):]
    parts = relative.rsplit('/', 1)
    if len(parts) != 2:
        return None
    source_subpath, shard_filename = parts
    match = SHARD_RE.match(shard_filename)
    if not match:
        return None
    traj_type, split = match.group(1), match.group(2)
    flat_subpath = source_subpath.replace('/', '-')
    return traj_type, split, flat_subpath, shard_filename


def group_keys(keys, source_prefix):
    """Group keys by (traj_type, split). Returns dict of (type, split) -> list of (key, flat_subpath)."""
    groups = defaultdict(list)
    skipped = []
    for key in keys:
        classified = classify_key(key, source_prefix)
        if classified is None:
            skipped.append(key)
            continue
        traj_type, split, flat_subpath, _ = classified
        groups[(traj_type, split)].append((key, flat_subpath))
    return groups, skipped


def download_keys(s3, bucket, keys, staging_dir):
    """Download all keys into staging_dir, preserving key paths under it.

    Returns the local path corresponding to each key.
    """
    paths = {}
    with tqdm(total=len(keys), desc="Downloading", unit="file") as pbar:
        for key in keys:
            local_path = os.path.join(staging_dir, key)
            os.makedirs(os.path.dirname(local_path), exist_ok=True)
            if not os.path.exists(local_path):
                s3.download_file(bucket, key, local_path)
            paths[key] = local_path
            pbar.update(1)
    return paths


def iter_tarball_members(path):
    """Yield (member_name, content_bytes) for each file member of a tar.zst archive."""
    dctx = zstd.ZstdDecompressor()
    with open(path, 'rb') as raw:
        with dctx.stream_reader(raw) as stream:
            with tarfile.open(fileobj=stream, mode='r|') as tar:
                for member in tar:
                    if not member.isfile():
                        continue
                    fileobj = tar.extractfile(member)
                    if fileobj is None:
                        continue
                    yield member.name, fileobj.read()


def count_fasta_frames_bases(content_str):
    """Return (frames, bases) for a FASTA string. Frames = header count; bases = total ACGT/N/- chars."""
    frames = 0
    bases = 0
    for line in content_str.splitlines():
        if line.startswith('>'):
            frames += 1
        else:
            bases += len(line.strip())
    return frames, bases


def build_pooled(s3, bucket, source_prefix, output_dir, shard_size,
                 staging_dir, keep_staging, shuffle):
    """Run the full rollup pipeline. Returns the metadata dict."""
    print(f"Listing s3://{bucket}/{source_prefix}/ ...")
    all_keys = list_source_keys(s3, bucket, source_prefix)
    print(f"Found {len(all_keys)} candidate tar.zst objects")

    groups, skipped = group_keys(all_keys, source_prefix)
    if skipped:
        print(f"Skipping {len(skipped)} keys that don't match shard pattern:")
        for k in skipped[:5]:
            print(f"  {k}")
        if len(skipped) > 5:
            print(f"  ... and {len(skipped) - 5} more")

    if not groups:
        sys.exit("No matching tar.zst shards found under source prefix")

    print("Source key counts per (type, split):")
    for (typ, split), entries in sorted(groups.items()):
        print(f"  {typ}-{split}: {len(entries)} source shards")

    keys_to_download = sorted({key for entries in groups.values() for key, _ in entries})

    os.makedirs(staging_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    print(f"Staging into {staging_dir}")
    download_keys(s3, bucket, keys_to_download, staging_dir)

    metadata = {
        'source_bucket': bucket,
        'source_prefix': source_prefix.rstrip('/'),
        'generated_at': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        'git_commit': get_git_commit(),
        'shard_size': shard_size,
        'shuffle_seed': 42 if shuffle else None,
        'splits': {},
        'source_keys': sorted(
            key[len(source_prefix.rstrip('/') + '/'):]
            for key in keys_to_download
        ),
    }

    for (traj_type, split), entries in sorted(groups.items()):
        prefix_name = f"{traj_type}-{split}"
        print(f"\nProcessing {prefix_name} ({len(entries)} source shards) ...")
        stats = {'trajectories': 0, 'frames': 0, 'bases': 0,
                 'shards': 0, 'num_source_keys': len(entries)}
        sorted_entries = sorted(entries, key=lambda e: e[0])

        writer = ShardWriter(output_dir, prefix_name, shard_size=shard_size,
                             shuffle=shuffle)
        for key, flat_subpath in tqdm(sorted_entries, desc=prefix_name, unit="src"):
            staged_path = os.path.join(staging_dir, key)
            for member_name, content_bytes in iter_tarball_members(staged_path):
                content_str = content_bytes.decode('utf-8')
                new_name = f"{flat_subpath}-{member_name}"
                writer.add(new_name, content_str)
                stats['trajectories'] += 1
                frames, bases = count_fasta_frames_bases(content_str)
                stats['frames'] += frames
                stats['bases'] += bases
        num_shards, num_files, _ = writer.close()
        stats['shards'] = num_shards
        assert num_files == stats['trajectories'], (
            f"writer counted {num_files} files but we tracked {stats['trajectories']}"
        )

        metadata['splits'][prefix_name] = stats
        print(f"  -> {stats['trajectories']} trajectories, {stats['frames']} frames, "
              f"{stats['bases']} bases, {stats['shards']} shards")

    metadata_path = os.path.join(output_dir, 'metadata.json')
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    print(f"\nWrote {metadata_path}")

    if not keep_staging:
        shutil.rmtree(staging_dir)
        print(f"Removed staging dir {staging_dir}")
    else:
        print(f"Kept staging dir {staging_dir}")

    return metadata


def dry_run(s3, bucket, source_prefix):
    """List + classify source keys without downloading or writing anything."""
    print(f"Listing s3://{bucket}/{source_prefix}/ ...")
    all_keys = list_source_keys(s3, bucket, source_prefix)
    print(f"Found {len(all_keys)} candidate tar.zst objects (pooled/ excluded)")

    groups, skipped = group_keys(all_keys, source_prefix)

    print("\nSource shards per (type, split):")
    if not groups:
        print("  (none)")
    for (typ, split), entries in sorted(groups.items()):
        unique_subpaths = {flat for _, flat in entries}
        print(f"  {typ}-{split}: {len(entries)} shards across {len(unique_subpaths)} sources")

    if skipped:
        print(f"\nSkipped {len(skipped)} keys not matching shard pattern:")
        for k in skipped[:10]:
            print(f"  {k}")
        if len(skipped) > 10:
            print(f"  ... and {len(skipped) - 10} more")

    print(f"\nTotal matched source shards: {sum(len(e) for e in groups.values())}")


def print_summary(metadata):
    print("\n=== Summary ===")
    print(f"Source: s3://{metadata['source_bucket']}/{metadata['source_prefix']}/")
    print(f"Generated at: {metadata['generated_at']}")
    print(f"Git commit: {metadata['git_commit']}")
    total_trajs = sum(s['trajectories'] for s in metadata['splits'].values())
    total_shards = sum(s['shards'] for s in metadata['splits'].values())
    print(f"{'Split':<20} {'Trajectories':>14} {'Frames':>14} {'Bases':>16} {'Shards':>8}")
    for name, stats in sorted(metadata['splits'].items()):
        print(f"{name:<20} {stats['trajectories']:>14,} {stats['frames']:>14,} "
              f"{stats['bases']:>16,} {stats['shards']:>8,}")
    print(f"{'TOTAL':<20} {total_trajs:>14,} {'':>14} {'':>16} {total_shards:>8,}")


def main():
    parser = argparse.ArgumentParser(
        description="Repackage granular S3 tar.zst trajectory shards into pooled shards."
    )
    parser.add_argument(
        "--bucket", default=None,
        help="S3 bucket (defaults to $S3_BUCKET env var)"
    )
    parser.add_argument(
        "--source-prefix", required=True,
        help="S3 key prefix to walk (e.g. 'trellis-trajectories' or 'trajectories/viral')"
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Local directory for pooled shards + metadata.json"
    )
    parser.add_argument(
        "--shard-size", type=int, default=10000,
        help="Trajectories per pooled shard (default: 10000)"
    )
    parser.add_argument(
        "--staging-dir", default=None,
        help="Local directory to stage downloads (default: <output_dir>/.staging)"
    )
    parser.add_argument(
        "--keep-staging", action="store_true",
        help="Keep staged downloads after processing (default: remove)"
    )
    parser.add_argument(
        "--no-shuffle", action="store_true",
        help="Disable in-shard shuffle (default: shuffle with seed=42)"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="List + group source keys, print summary, and exit (no downloads or writes)"
    )
    args = parser.parse_args()

    bucket = args.bucket or os.environ.get('S3_BUCKET')
    if not bucket:
        sys.exit("Error: --bucket not provided and $S3_BUCKET is not set")

    for var in ('AWS_ACCESS_KEY_ID', 'AWS_SECRET_ACCESS_KEY'):
        if not os.environ.get(var):
            sys.exit(f"Error: Missing required environment variable: {var}")

    s3 = boto3.client('s3')

    if args.dry_run:
        dry_run(s3, bucket, args.source_prefix)
        return

    staging_dir = args.staging_dir or os.path.join(args.output_dir, '.staging')
    metadata = build_pooled(
        s3=s3,
        bucket=bucket,
        source_prefix=args.source_prefix,
        output_dir=args.output_dir,
        shard_size=args.shard_size,
        staging_dir=staging_dir,
        keep_staging=args.keep_staging,
        shuffle=not args.no_shuffle,
    )
    print_summary(metadata)


if __name__ == "__main__":
    main()
