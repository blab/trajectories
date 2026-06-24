"""Upload results directory to S3 bucket."""

import argparse
import os
import re
import sys
import boto3
from tqdm import tqdm


EXCLUDED_FILES = {'.DS_Store', '.snakemake_timestamp'}


def get_s3_prefix_for_analysis(analysis, base_prefix):
    """Get S3 prefix, nesting subtrees under their parent family.

    For example:
      - rdrp-paramyxoviridae-xs -> trajectories/rdrp-paramyxoviridae-xs/
      - rdrp-paramyxoviridae-xs_001 -> trajectories/rdrp-paramyxoviridae-xs/subtrees/rdrp-paramyxoviridae-xs_001/
      - cytb-xs -> trajectories/cytb-xs/
    """
    # Match rdrp-FAMILY-xs_XXX pattern (subtrees)
    match = re.match(r'^(rdrp-[a-z]+-xs)_(\d+)$', analysis)
    if match:
        parent = match.group(1)
        # Subtrees go under their parent
        return f"{base_prefix}/{parent}/subtrees/{analysis}"
    return f"{base_prefix}/{analysis}"


def count_files(directory):
    """Count .tar.zst shards in directory tree."""
    total = 0
    for root, dirs, files in os.walk(directory):
        total += sum(1 for f in files if f.endswith('.tar.zst'))
    return total


def delete_s3_prefix(s3, bucket, prefix, exclude_subdir=None):
    """Delete all objects under an S3 prefix, optionally excluding a subdirectory."""
    paginator = s3.get_paginator('list_objects_v2')
    deleted_count = 0
    exclude_prefix = f"{prefix}{exclude_subdir}/" if exclude_subdir else None
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        if 'Contents' in page:
            objects = []
            for obj in page['Contents']:
                # Skip objects in excluded subdirectory
                if exclude_prefix and obj['Key'].startswith(exclude_prefix):
                    continue
                objects.append({'Key': obj['Key']})
            if objects:
                s3.delete_objects(Bucket=bucket, Delete={'Objects': objects})
                deleted_count += len(objects)
    return deleted_count


def main():
    parser = argparse.ArgumentParser(
        description="Upload results directory to S3 bucket."
    )
    parser.add_argument(
        "--upload-dir", default="results",
        help="Local directory to upload"
    )
    parser.add_argument(
        "--prefix", default="trajectories",
        help="S3 key prefix"
    )
    parser.add_argument(
        "--filter", default=None,
        help="Regex pattern to filter which datasets to upload (e.g., 'rdrp-.*-xs$' for main RdRp only)"
    )
    parser.add_argument(
        "--analyses", nargs="*", default=None,
        help="Explicit list of analysis names to upload (takes precedence over --filter)"
    )
    parser.add_argument(
        "--branches-dir", default="data",
        help="Directory holding {analysis}/branches.tsv to upload alongside the shards "
             "(the tree edge-list + train/test split; trajectories->tree is not invertible)"
    )
    args = parser.parse_args()

    # Check required environment variables
    required_vars = ['S3_BUCKET', 'AWS_ACCESS_KEY_ID', 'AWS_SECRET_ACCESS_KEY']
    missing = [v for v in required_vars if not os.environ.get(v)]
    if missing:
        sys.exit(f"Error: Missing required environment variables: {', '.join(missing)}")

    bucket = os.environ['S3_BUCKET']

    # Initialize S3 client (uses AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY from env)
    s3 = boto3.client('s3')

    upload_dir = args.upload_dir

    # Compile filter pattern if provided
    filter_pattern = re.compile(args.filter) if args.filter else None

    # Find datasets to upload
    if args.analyses is not None:
        # Explicit list of analyses takes precedence
        datasets_to_upload = [a for a in args.analyses if os.path.isdir(os.path.join(upload_dir, a))]
        missing = [a for a in args.analyses if not os.path.isdir(os.path.join(upload_dir, a))]
        if missing:
            print(f"Warning: analysis directories not found: {', '.join(missing)}")
        print(f"Uploading {len(datasets_to_upload)} specified analysis dataset(s): {', '.join(sorted(datasets_to_upload))}")
    else:
        # Fall back to directory scan with optional filter
        datasets_to_upload = []
        for item in os.listdir(upload_dir):
            item_path = os.path.join(upload_dir, item)
            if os.path.isdir(item_path):
                if filter_pattern and not filter_pattern.match(item):
                    continue
                datasets_to_upload.append(item)
        if filter_pattern:
            print(f"Filter '{args.filter}' matched {len(datasets_to_upload)} dataset(s): {', '.join(sorted(datasets_to_upload))}")

    # Delete existing S3 objects for each dataset directory before uploading
    # This prevents stale files from accumulating when dataset tips change
    for item in datasets_to_upload:
        s3_prefix = get_s3_prefix_for_analysis(item, args.prefix) + "/"
        # For main datasets (not subtrees), exclude the subtrees/ subdirectory from deletion
        is_subtree = re.match(r'^rdrp-[a-z]+-xs_\d+$', item)
        exclude_subdir = None if is_subtree else "subtrees"
        print(f"Deleting existing files at s3://{bucket}/{s3_prefix}" +
              (f" (excluding {exclude_subdir}/)" if exclude_subdir else ""))
        deleted = delete_s3_prefix(s3, bucket, s3_prefix, exclude_subdir=exclude_subdir)
        if deleted:
            print(f"  Deleted {deleted} existing files")

    # Count files for progress bar (only in datasets we're uploading)
    total_files = 0
    for dataset in datasets_to_upload:
        total_files += count_files(os.path.join(upload_dir, dataset))
        if os.path.isfile(os.path.join(args.branches_dir, dataset, "branches.tsv")):
            total_files += 1
    # Also count root-level files (like summary.json) unless using --filter
    if not filter_pattern:
        for f in os.listdir(upload_dir):
            if os.path.isfile(os.path.join(upload_dir, f)) and f not in EXCLUDED_FILES:
                total_files += 1

    print(f"Uploading {total_files} files to s3://{bucket}/{args.prefix}/")

    # Upload files from filtered datasets
    with tqdm(total=total_files, desc="Uploading") as pbar:
        # Upload root-level files (like summary.json) unless using --filter
        if not filter_pattern:
            for filename in os.listdir(upload_dir):
                local_path = os.path.join(upload_dir, filename)
                if os.path.isfile(local_path) and filename not in EXCLUDED_FILES:
                    s3_key = f"{args.prefix}/{filename}"
                    s3.upload_file(local_path, bucket, s3_key)
                    pbar.update(1)

        # Upload dataset directories
        for dataset in datasets_to_upload:
            dataset_path = os.path.join(upload_dir, dataset)
            for root, dirs, files in os.walk(dataset_path):
                for filename in files:
                    if not filename.endswith('.tar.zst'):
                        continue
                    local_path = os.path.join(root, filename)
                    # Get path relative to the dataset directory
                    relative_to_dataset = os.path.relpath(local_path, dataset_path)
                    # Get the proper S3 prefix for this analysis (handles subtree nesting)
                    s3_analysis_prefix = get_s3_prefix_for_analysis(dataset, args.prefix)
                    s3_key = f"{s3_analysis_prefix}/{relative_to_dataset}"
                    s3.upload_file(local_path, bucket, s3_key)
                    pbar.update(1)

            # Upload branches.tsv (parent/child/hamming + train_test) from the data
            # dir into the same prefix. The trajectories are a lossy, truncated sample
            # of this tree (trajectories->tree is not invertible), so shipping the
            # edge-list lets consumers recover the full tree + the exact train/test split.
            branches_path = os.path.join(args.branches_dir, dataset, "branches.tsv")
            if os.path.isfile(branches_path):
                s3_analysis_prefix = get_s3_prefix_for_analysis(dataset, args.prefix)
                s3.upload_file(branches_path, bucket, f"{s3_analysis_prefix}/branches.tsv")
                pbar.update(1)

    print(f"Done! Uploaded to s3://{bucket}/{args.prefix}/")


if __name__ == "__main__":
    main()
