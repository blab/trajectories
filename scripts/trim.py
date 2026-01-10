# code from ChatGPT
import argparse
from Bio import SeqIO

def trim_alignment(input_file, output_file, begin, end):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if begin is not None and end is not None:
                # Trim to specified positions
                trimmed_sequence = record.seq[begin-1:end]  # Convert to 0-based indexing
            else:
                # No trimming - keep full sequence
                trimmed_sequence = record.seq
            record.seq = trimmed_sequence
            SeqIO.write(record, outfile, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Trim an alignment to specified positions.")
    parser.add_argument("--input-alignment", type=str, required=True, help="Path to the input alignment in FASTA format.")
    parser.add_argument("--output-alignment", type=str, default="data/trimmed.fasta", help="Path to save the trimmed alignment.")
    parser.add_argument("--begin", type=int, default=None, help="Start position (1-based, inclusive). If not provided, no trimming is done.")
    parser.add_argument("--end", type=int, default=None, help="End position (1-based, inclusive). If not provided, no trimming is done.")

    args = parser.parse_args()

    if args.begin is not None and args.end is not None:
        print(f"Trimming alignment from positions {args.begin} to {args.end}...")
    else:
        print("No trimming specified, copying full sequences...")
    trim_alignment(args.input_alignment, args.output_alignment, args.begin, args.end)
    print(f"Output saved to {args.output_alignment}")

if __name__ == "__main__":
    main()
