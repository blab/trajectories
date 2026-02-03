"""
ShardWriter: Buffer trajectories in memory and write directly to tar.zst shards.

Eliminates filesystem overhead from writing millions of individual .fasta files
by streaming directly to compressed tar archives.
"""

import io
import os
import random
import tarfile
import zstandard as zstd


class ShardWriter:
    """
    Buffers trajectories and writes directly to tar.zst shards.

    Usage:
        with ShardWriter("output/", "forwards-train", shard_size=10000) as writer:
            for tip in tips:
                content = build_trajectory(tip)
                writer.add(f"{tip}.fasta", content)

        # After context exits, writer.stats contains (shards, files, bytes)
    """

    def __init__(self, output_dir, prefix, shard_size=10000, shuffle=True, seed=42):
        """
        Initialize ShardWriter.

        Args:
            output_dir: Directory to write shard files to
            prefix: Prefix for shard filenames (e.g., "forwards-train")
            shard_size: Number of files per shard
            shuffle: Whether to shuffle files within each shard before writing
            seed: Random seed for shuffling
        """
        self.output_dir = output_dir
        self.prefix = prefix
        self.shard_size = shard_size
        self.shuffle = shuffle
        self.seed = seed

        self.buffer = []  # List of (filename, content) tuples
        self.shard_idx = 0
        self.total_files = 0
        self.total_bytes = 0

        self.compressor = zstd.ZstdCompressor()
        self.rng = random.Random(seed)

        # Create output directory
        os.makedirs(output_dir, exist_ok=True)

    def add(self, filename, content):
        """
        Add a file to the buffer.

        Args:
            filename: Name for the file within the tar archive
            content: String content of the file
        """
        self.buffer.append((filename, content))
        if len(self.buffer) >= self.shard_size:
            self._flush_shard()

    def _flush_shard(self):
        """Shuffle buffer (if enabled), write to tar.zst, clear buffer."""
        if not self.buffer:
            return

        if self.shuffle:
            self.rng.shuffle(self.buffer)

        # Create tar archive in memory
        tar_buffer = io.BytesIO()
        with tarfile.open(fileobj=tar_buffer, mode='w') as tar:
            for filename, content in self.buffer:
                # Create TarInfo for the file
                data = content.encode('utf-8')
                info = tarfile.TarInfo(name=filename)
                info.size = len(data)
                tar.addfile(info, io.BytesIO(data))

        # Compress with zstd and write
        tar_data = tar_buffer.getvalue()
        compressed = self.compressor.compress(tar_data)

        # Write shard file
        shard_path = os.path.join(
            self.output_dir,
            f"{self.prefix}-{self.shard_idx:03d}.tar.zst"
        )
        with open(shard_path, 'wb') as f:
            f.write(compressed)

        # Update stats
        self.total_files += len(self.buffer)
        self.total_bytes += len(compressed)
        self.shard_idx += 1

        # Clear buffer
        self.buffer = []

    def close(self):
        """
        Flush remaining buffer and return stats.

        Returns:
            Tuple of (num_shards, num_files, total_bytes)
        """
        self._flush_shard()
        return (self.shard_idx, self.total_files, self.total_bytes)

    @property
    def stats(self):
        """Return current stats (may change if not closed)."""
        return (self.shard_idx, self.total_files, self.total_bytes)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False
