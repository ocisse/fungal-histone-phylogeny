#!/bin/bash

# Exit on any error
set -e

# === CONFIGURE THESE PATHS ===
BASE_DIR=~/PCP_genomics/KINETO/data/processed/phylogeny/histones/motifs/clustering
IUPRED_DIR=~/PCP_genomics/KINETO/src/tools/iupred3  # directory containing iupred3.py and iupred3_lib.py

# === FIXED PATHS ===
INPUT_DIR="$BASE_DIR/cdhit_cluster_fastas"
OUTPUT_DIR="$BASE_DIR/iupred3_results"
LOG_FILE="$BASE_DIR/iupred3_global_scores.log"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Create or clear log file
echo "# Cluster : Global_Disorder_Fraction" > "$LOG_FILE"

# Loop over each cluster FASTA file
for fasta in "$INPUT_DIR"/*.fasta; do
  base=$(basename "$fasta" .fasta)
  out_file="$OUTPUT_DIR/${base}_iupred.txt"

  echo "Running IUPred3 on $base..."

  # Run IUPred3 using full path, and set PYTHONPATH to find the iupred3_lib module
  PYTHONPATH="$IUPRED_DIR" python "$IUPRED_DIR/iupred3.py" "$fasta" long > "$out_file"

  # Compute global disorder score (average fraction of disordered residues across all sequences)
  awk '
    BEGIN { FS="\t" }
    /^[^#]/ {
      id = $1
      score = $3
      total[id]++
      if (score >= 0.5) disordered[id]++
    }
    END {
      for (id in total) {
        fraction = disordered[id] / total[id]
        sum += fraction
        count++
      }
      if (count > 0)
        printf "%s : %.4f\n", "'$base'", sum / count
      else
        printf "%s : NA\n", "'$base'"
    }
  ' "$out_file" >> "$LOG_FILE"
done

echo "Done running IUPred3 on all clusters."
echo "Global disorder scores saved in: $LOG_FILE"

