#!/usr/bin/env python3

import os
from collections import defaultdict

# File paths (adjust if needed)
REAL_SCORE_FILE = "/PCP_genomics/KINETO/data/processed/phylogeny/histones/motifs/clustering/glam2_real_scores.log"
SHUFFLED_SCORE_FILE = "/PCP_genomics/KINETO/data/processed/phylogeny/histones/motifs/clustering/glam2_shuffled_scores.log"
OUTPUT_FILE = "/PCP_genomics/KINETO/data/processed/phylogeny/histones/motifs/clustering/glam2_empirical_pvalues.tsv"

def read_scores(filepath, is_shuffled=False):
    scores = defaultdict(list) if is_shuffled else {}
    with open(filepath, "r") as f:
        prev_cluster = None
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            if ":" in line:
                cluster, score = line.split(":")
                cluster = cluster.strip().replace(".fasta.glam2_out", "")
                try:
                    score_val = float(score.strip())
                    if is_shuffled:
                        base = cluster.rsplit("_shuffled_", 1)[0]
                        scores[base].append(score_val)
                    else:
                        scores[cluster] = score_val
                    prev_cluster = cluster
                except ValueError:
                    continue
            elif prev_cluster and not is_shuffled:
                # If this line is a second line with a score (without a label), ignore it
                continue
    return scores

# Read real and shuffled scores
real_scores = read_scores(REAL_SCORE_FILE, is_shuffled=False)
shuffled_scores = read_scores(SHUFFLED_SCORE_FILE, is_shuffled=True)

# Write output
with open(OUTPUT_FILE, "w") as out:
    out.write("Cluster\tRealScore\tShuffledMean\tEmpiricalPvalue\n")
    for cluster, real_score in real_scores.items():
        shuffles = shuffled_scores.get(cluster, [])
        if not shuffles:
            continue
        count = sum(1 for s in shuffles if s >= real_score)
        pval = (count + 1) / (len(shuffles) + 1)  # smoothed empirical p-value
        mean_shuffle = sum(shuffles) / len(shuffles)
        out.write(f"{cluster}\t{real_score:.3f}\t{mean_shuffle:.3f}\t{pval:.4g}\n")

print(f"Done. Empirical p-values saved to: {OUTPUT_FILE}")

