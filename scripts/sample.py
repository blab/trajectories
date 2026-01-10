#!/usr/bin/env python3
"""
Sample a subset of sequences from a FASTA file for faster testing.
"""

import argparse
import random
from Bio import SeqIO

def sample_sequences(input_file, output_file, sample_fraction=0.1, include_nodes=False):
    """Sample a fraction of sequences from input FASTA and save to output FASTA."""
    
    # Load all sequences
    sequences = list(SeqIO.parse(input_file, "fasta"))
    
    # Filter sequences based on type (tips vs nodes)
    if include_nodes:
        # Include all sequences (tips and internal nodes)
        filtered_sequences = sequences
        sequence_type = "all sequences"
    else:
        # Only include tip sequences (exclude internal nodes that start with NODE_)
        filtered_sequences = [seq for seq in sequences if not seq.id.startswith('NODE_')]
        sequence_type = "tip sequences only"
    
    total_sequences = len(filtered_sequences)
    
    if total_sequences == 0:
        raise ValueError(f"No sequences found after filtering for {sequence_type}")
    
    # Calculate number of samples
    n_samples = int(total_sequences * sample_fraction)
    n_samples = max(1, n_samples)  # Ensure at least 1 sequence
    
    # Set random seed for reproducible sampling
    random.seed(42)
    
    # Sample sequences
    sampled_sequences = random.sample(filtered_sequences, n_samples)
    
    # Write sampled sequences
    SeqIO.write(sampled_sequences, output_file, "fasta")
    
    print(f"Filtered to {total_sequences} {sequence_type} from {len(sequences)} total sequences")
    print(f"Sampled {n_samples} out of {total_sequences} sequences ({sample_fraction:.1%})")
    print(f"Sampled sequences saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sample a subset of sequences from a FASTA file")
    parser.add_argument("--input", type=str, required=True, help="Path to the input FASTA file")
    parser.add_argument("--output", type=str, required=True, help="Path to save the sampled FASTA file")
    parser.add_argument("--fraction", type=float, default=0.1, help="Fraction of sequences to sample (0.0-1.0, default: 0.1)")
    parser.add_argument("--include-nodes", action="store_true", help="Include internal nodes (NODE_*) in sampling, otherwise sample only tips")
    
    args = parser.parse_args()
    
    if not 0.0 < args.fraction <= 1.0:
        raise ValueError("Sample fraction must be between 0.0 and 1.0")
    
    sample_sequences(args.input, args.output, args.fraction, args.include_nodes)