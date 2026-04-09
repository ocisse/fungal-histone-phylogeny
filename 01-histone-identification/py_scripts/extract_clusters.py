#!/usr/bin/env python3

import argparse
import os
import re
from collections import defaultdict
from Bio import SeqIO

def sanitize_filename(name):
    return re.sub(r"[^\w\-_.]", "_", name)

def load_clusters(cluster_file):
    clusters = defaultdict(list)
    with open(cluster_file, "r") as f:
        for line in f:
            rep, member = line.strip().split('\t')
            clusters[rep].append(member)
    return clusters

def main():
    parser = argparse.ArgumentParser(description="Extract per-cluster FASTA files from MMseqs2 clustering.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file used in clustering")
    parser.add_argument("-c", "--clusters", required=True, help="MMseqs2 result_cluster.tsv file")
    parser.add_argument("-o", "--outdir", default="clusters_out", help="Output directory (default: clusters_out)")
    parser.add_argument("--min-size", type=int, default=1, help="Minimum cluster size to output (default: 1)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("[INFO] Reading input FASTA...")
    seqs = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))

    print("[INFO] Parsing cluster assignments...")
    clusters = load_clusters(args.clusters)

    print(f"[INFO] Writing clusters with at least {args.min_size} members...")
    written = 0
    for rep, members in clusters.items():
        if len(members) < args.min_size:
            continue
        safe_rep = sanitize_filename(rep)
        outfile = os.path.join(args.outdir, f"cluster_{safe_rep}.fasta")
        with open(outfile, "w") as out_f:
            for m in members:
                if m in seqs:
                    SeqIO.write(seqs[m], out_f, "fasta")
                else:
                    print(f"[WARNING] Sequence '{m}' not found in input FASTA.")
        written += 1

    print(f"[DONE] {written} cluster FASTA files written to '{args.outdir}'.")

if __name__ == "__main__":
    main()

