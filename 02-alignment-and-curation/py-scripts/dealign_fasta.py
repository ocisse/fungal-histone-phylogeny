#!/usr/bin/env python3

from Bio import SeqIO
import sys

def dealign_fasta(input_fasta, output_fasta):
    with open(input_fasta, "r") as in_handle, open(output_fasta, "w") as out_handle:
        for record in SeqIO.parse(in_handle, "fasta"):
            # Remove gaps from the sequence
            record.seq = record.seq.ungap("-")
            SeqIO.write(record, out_handle, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python dealign_fasta.py <input.fasta> <output.fasta>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]

    dealign_fasta(input_fasta, output_fasta)
    print(f"De-aligned FASTA written to: {output_fasta}")

