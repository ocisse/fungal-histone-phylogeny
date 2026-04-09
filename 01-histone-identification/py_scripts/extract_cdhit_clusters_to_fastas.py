#!/usr/bin/env python3

import os
import sys
import re
import argparse
from collections import defaultdict
from Bio import SeqIO

def parse_clstr(clstr_file):
    clusters = defaultdict(list)
    current_cluster = None
    pattern = re.compile(r'>Spec\d+\|[^\s\.]+')

    with open(clstr_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                current_cluster = line.split()[1]
            else:
                match = pattern.search(line)
                if match:
                    seq_id = match.group(0)[1:]  # remove '>'
                    clusters[current_cluster].append(seq_id)
    return clusters

def write_cluster_fastas(clusters, fasta_file, output_dir, min_size=1):
    os.makedirs(output_dir, exist_ok=True)

    fasta_records = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    fasta_ids = list(fasta_records.keys())
    count = 0

    for cluster_id, seq_ids in clusters.items():
        if len(seq_ids) < min_size:
            continue  # skip small clusters

        output_path = os.path.join(output_dir, f"cluster_{cluster_id}.fasta")
        with open(output_path, "w") as out_f:
            for seq_id in seq_ids:
                # Try exact match first
                if seq_id in fasta_records:
                    SeqIO.write(fasta_records[seq_id], out_f, "fasta")
                else:
                    # fallback: find fasta ID starting with seq_id
                    matched_id = next((fid for fid in fasta_ids if fid.startswith(seq_id)), None)
                    if matched_id:
                        SeqIO.write(fasta_records[matched_id], out_f, "fasta")
                    else:
                        print(f"⚠️  Warning: sequence ID '{seq_id}' not found in FASTA", file=sys.stderr)
        count += 1

    print(f"✅ Extracted {count} clusters with ≥{min_size} sequences to: {output_dir}")

def main():
    parser = argparse.ArgumentParser(description="Extract CD-HIT clusters to individual FASTA files.")
    parser.add_argument("clstr_file", help="CD-HIT .clstr file")
    parser.add_argument("fasta_file", help="Original FASTA file")
    parser.add_argument("output_dir", help="Output directory for cluster FASTAs")
    parser.add_argument("--min-size", type=int, default=1, help="Minimum cluster size to include (default: 1)")

    args = parser.parse_args()

    clusters = parse_clstr(args.clstr_file)
    write_cluster_fastas(clusters, args.fasta_file, args.output_dir, min_size=args.min_size)

if __name__ == "__main__":
    main()

