"""Package trajectory FASTA files into sharded tar.zst archives."""

import argparse
import io
import os
import random
import tarfile
import zstandard as zstd
from tqdm import tqdm


def get_fasta_files(input_dir):
    """Get list of .fasta files in directory."""
    files = []
    for filename in os.listdir(input_dir):
        if filename.endswith('.fasta'):
            files.append(filename)
    return sorted(files)


def create_shard(files, input_dir, output_path):
    """Create a tar.zst shard from a list of FASTA files."""
    # Create tar archive in memory
    tar_buffer = io.BytesIO()
    with tarfile.open(fileobj=tar_buffer, mode='w') as tar:
        for filename in files:
            filepath = os.path.join(input_dir, filename)
            tar.add(filepath, arcname=filename)

    # Compress with zstd and write
    tar_data = tar_buffer.getvalue()
    cctx = zstd.ZstdCompressor()
    compressed = cctx.compress(tar_data)

    with open(output_path, 'wb') as f:
        f.write(compressed)

    return len(tar_data), len(compressed)


def main():
    parser = argparse.ArgumentParser(
        description="Package trajectory FASTA files into sharded tar.zst archives."
    )
    parser.add_argument(
        "--input-dir", required=True,
        help="Directory containing .fasta files"
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Output directory for tar.zst shards"
    )
    parser.add_argument(
        "--shard-size", type=int, default=10000,
        help="Number of trajectories per shard (default: 10000)"
    )
    parser.add_argument(
        "--shuffle", action="store_true",
        help="Shuffle files before sharding"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for shuffling (default: 42)"
    )
    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Get list of FASTA files
    files = get_fasta_files(args.input_dir)
    if not files:
        print(f"No .fasta files found in {args.input_dir}")
        return

    print(f"Found {len(files)} FASTA files")

    # Shuffle if requested
    if args.shuffle:
        random.seed(args.seed)
        random.shuffle(files)
        print(f"Shuffled files (seed={args.seed})")

    # Group into shards
    shards = []
    for i in range(0, len(files), args.shard_size):
        shards.append(files[i:i + args.shard_size])

    print(f"Creating {len(shards)} shard(s) with up to {args.shard_size} files each")

    # Create each shard
    total_uncompressed = 0
    total_compressed = 0
    for shard_idx, shard_files in enumerate(tqdm(shards, desc="Creating shards")):
        output_path = os.path.join(
            args.output_dir,
            f"trajectories-{shard_idx:03d}.tar.zst"
        )
        uncompressed, compressed = create_shard(shard_files, args.input_dir, output_path)
        total_uncompressed += uncompressed
        total_compressed += compressed

    # Report compression stats
    ratio = total_uncompressed / total_compressed if total_compressed > 0 else 0
    print(f"Done! Created {len(shards)} shard(s) in {args.output_dir}")
    print(f"Total size: {total_compressed / 1024 / 1024:.1f} MB "
          f"(compression ratio: {ratio:.1f}x)")


if __name__ == "__main__":
    main()
