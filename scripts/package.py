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


def process_split(input_dir, output_dir, shard_size, shuffle, seed, split_name=None):
    """
    Process FASTA files from a directory into sharded tar.zst archives.

    Args:
        input_dir: Directory containing .fasta files
        output_dir: Output directory for tar.zst shards
        shard_size: Number of trajectories per shard
        shuffle: Whether to shuffle files before sharding
        seed: Random seed for shuffling
        split_name: Optional split name (e.g., 'train', 'test') for shard naming

    Returns:
        Tuple of (num_shards, total_uncompressed, total_compressed)
    """
    files = get_fasta_files(input_dir)
    if not files:
        return 0, 0, 0

    print(f"Found {len(files)} {split_name or ''} FASTA files".strip())

    if shuffle:
        random.seed(seed)
        random.shuffle(files)

    # Group into shards
    shards = []
    for i in range(0, len(files), shard_size):
        shards.append(files[i:i + shard_size])

    # Create each shard
    total_uncompressed = 0
    total_compressed = 0
    desc = f"Creating {split_name} shards" if split_name else "Creating shards"
    for shard_idx, shard_files in enumerate(tqdm(shards, desc=desc)):
        if split_name:
            filename = f"forwards-{split_name}-{shard_idx:03d}.tar.zst"
        else:
            filename = f"trajectories-{shard_idx:03d}.tar.zst"
        output_path = os.path.join(output_dir, filename)
        uncompressed, compressed = create_shard(shard_files, input_dir, output_path)
        total_uncompressed += uncompressed
        total_compressed += compressed

    return len(shards), total_uncompressed, total_compressed


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

    if args.shuffle:
        print(f"Shuffling enabled (seed={args.seed})")

    # Check for train/test subdirectories
    train_dir = os.path.join(args.input_dir, 'forwards-train')
    test_dir = os.path.join(args.input_dir, 'forwards-test')

    total_shards = 0
    total_uncompressed = 0
    total_compressed = 0

    if os.path.isdir(train_dir) and os.path.isdir(test_dir):
        # Process train and test separately
        for split in ['train', 'test']:
            split_dir = os.path.join(args.input_dir, f'forwards-{split}')
            num_shards, uncompressed, compressed = process_split(
                split_dir, args.output_dir, args.shard_size,
                args.shuffle, args.seed, split_name=split
            )
            total_shards += num_shards
            total_uncompressed += uncompressed
            total_compressed += compressed
    else:
        # Fall back to original behavior (backward compatibility)
        num_shards, uncompressed, compressed = process_split(
            args.input_dir, args.output_dir, args.shard_size,
            args.shuffle, args.seed
        )
        total_shards = num_shards
        total_uncompressed = uncompressed
        total_compressed = compressed

    if total_shards == 0:
        print(f"No .fasta files found in {args.input_dir}")
        return

    # Report compression stats
    ratio = total_uncompressed / total_compressed if total_compressed > 0 else 0
    print(f"Done! Created {total_shards} shard(s) in {args.output_dir}")
    print(f"Total size: {total_compressed / 1024 / 1024:.1f} MB "
          f"(compression ratio: {ratio:.1f}x)")


if __name__ == "__main__":
    main()
