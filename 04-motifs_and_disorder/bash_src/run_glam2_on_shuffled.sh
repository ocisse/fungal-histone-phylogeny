#!/bin/bash

# Exit on error
set -e

# Directories
BASE_DIR=~/PCP_genomics/KINETO/data/processed/phylogeny/histones/motifs/clustering
INPUT_DIR="$BASE_DIR/shuffled_cdhit_cluster_fastas"
OUTPUT_DIR="$BASE_DIR/glam2_shuffled_results"
LOG_FILE="$BASE_DIR/glam2_shuffled_scores.log"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"
echo "# File : Motif Score" > "$LOG_FILE"

# Run GLAM2 on each shuffled FASTA file
for fasta in "$INPUT_DIR"/*.fasta; do
  base=$(basename "$fasta" .fasta)
  out_dir="$OUTPUT_DIR/$base"
  mkdir -p "$out_dir"

  echo "Running GLAM2 on $base..."

  # Run glam2 (adjust parameters if needed)
  glam2 p -O "$out_dir" -n 1 "$fasta" > /dev/null

  # Extract motif score from output
  score=$(grep '^Score' "$out_dir/glam2.txt" | awk '{print $2}')

  # Log result
  echo "$base : $score" >> "$LOG_FILE"
done

echo "All GLAM2 runs complete."
echo "Scores saved in: $LOG_FILE"

