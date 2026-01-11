"""Upload results directory to S3 bucket."""

import argparse
import os
import sys
import boto3
from tqdm import tqdm


EXCLUDED_FILES = {'.DS_Store', '.snakemake_timestamp'}


def count_files(directory):
    """Count total files in directory tree (excluding system files)."""
    total = 0
    for root, dirs, files in os.walk(directory):
        total += sum(1 for f in files if f not in EXCLUDED_FILES)
    return total


def main():
    parser = argparse.ArgumentParser(
        description="Upload results directory to S3 bucket."
    )
    parser.add_argument(
        "--results-dir", default="results",
        help="Local results directory to upload"
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

    # Count files for progress bar
    total_files = count_files(args.results_dir)
    print(f"Uploading {total_files} files to s3://{bucket}/{args.prefix}/")

    # Walk results directory and upload all files
    with tqdm(total=total_files, desc="Uploading") as pbar:
        for root, dirs, files in os.walk(args.results_dir):
            for filename in files:
                if filename in EXCLUDED_FILES:
                    continue
                local_path = os.path.join(root, filename)
                # Convert local path to S3 key
                relative_path = os.path.relpath(local_path, args.results_dir)
                s3_key = f"{args.prefix}/{relative_path}"
                s3.upload_file(local_path, bucket, s3_key)
                pbar.update(1)

    print(f"Done! Uploaded to s3://{bucket}/{args.prefix}/")


if __name__ == "__main__":
    main()
