"""Upload export directory to S3 bucket."""

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
    """Count total files in directory tree (excluding system files)."""
    total = 0
    for root, dirs, files in os.walk(directory):
        total += sum(1 for f in files if f not in EXCLUDED_FILES)
    return total


def delete_s3_prefix(s3, bucket, prefix):
    """Delete all objects under an S3 prefix."""
    paginator = s3.get_paginator('list_objects_v2')
    deleted_count = 0
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        if 'Contents' in page:
            objects = [{'Key': obj['Key']} for obj in page['Contents']]
            s3.delete_objects(Bucket=bucket, Delete={'Objects': objects})
            deleted_count += len(objects)
    return deleted_count


def main():
    parser = argparse.ArgumentParser(
        description="Upload export directory to S3 bucket."
    )
    parser.add_argument(
        "--export-dir", default="export",
        help="Local export directory to upload"
    )
    parser.add_argument(
        "--prefix", default="trajectories",
        help="S3 key prefix"
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

    # Delete existing S3 objects for each dataset directory before uploading
    # This prevents stale files from accumulating when dataset tips change
    for item in os.listdir(args.export_dir):
        item_path = os.path.join(args.export_dir, item)
        if os.path.isdir(item_path):
            s3_prefix = get_s3_prefix_for_analysis(item, args.prefix) + "/"
            print(f"Deleting existing files at s3://{bucket}/{s3_prefix}")
            deleted = delete_s3_prefix(s3, bucket, s3_prefix)
            if deleted:
                print(f"  Deleted {deleted} existing files")

    # Count files for progress bar
    total_files = count_files(args.export_dir)
    print(f"Uploading {total_files} files to s3://{bucket}/{args.prefix}/")

    # Walk results directory and upload all files
    with tqdm(total=total_files, desc="Uploading") as pbar:
        for root, dirs, files in os.walk(args.export_dir):
            for filename in files:
                if filename in EXCLUDED_FILES:
                    continue
                local_path = os.path.join(root, filename)
                # Convert local path to S3 key with proper nesting
                relative_path = os.path.relpath(local_path, args.export_dir)
                # Get the analysis name (first path component)
                path_parts = relative_path.split(os.sep)

                # Handle root-level files (like summary.json)
                if len(path_parts) == 1:
                    s3_key = f"{args.prefix}/{filename}"
                else:
                    analysis = path_parts[0]
                    rest_of_path = os.sep.join(path_parts[1:])
                    # Get the proper S3 prefix for this analysis (handles subtree nesting)
                    s3_analysis_prefix = get_s3_prefix_for_analysis(analysis, args.prefix)
                    s3_key = f"{s3_analysis_prefix}/{rest_of_path}"

                s3.upload_file(local_path, bucket, s3_key)
                pbar.update(1)

    print(f"Done! Uploaded to s3://{bucket}/{args.prefix}/")


if __name__ == "__main__":
    main()
