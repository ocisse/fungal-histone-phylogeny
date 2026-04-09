#!/bin/bash

# Exit on any error
set -e

# Define directories
BASE_DIR=~/PCP_genomics/KINETO/data/processed/phylogeny/histones/motifs/clustering
INPUT_DIR="$BASE_DIR/cdhit_cluster_fastas"
OUTPUT_DIR="$BASE_DIR/glam2_real_results"
LOG_FILE="$BASE_DIR/glam2_real_scores.log"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Create or clear log file
echo "# File : Motif Score" > "$LOG_FILE"

# Loop over each real (clustered) FASTA file
for fasta in "$INPUT_DIR"/*.fasta; do
  base=$(basename "$fasta" .fasta)
  out_dir="$OUTPUT_DIR/$base"
  mkdir -p "$out_dir"

  echo "Running GLAM2 on $base..."

  # Run GLAM2 (adjust parameters as needed)
  glam2 -O "$out_dir" -n 1 "$fasta" > /dev/null

  # Extract motif score from glam2.txt
  score=$(grep '^Score' "$out_dir/glam2.txt" | awk '{print $2}')

  # Log result
  echo "$base : $score" >> "$LOG_FILE"
done

echo "Done running GLAM2 on all real sequences."
echo "Scores saved in: $LOG_FILE"

