#!/bin/bash

# Exit on any error
set -e

# Define GLAM2 output directory (adjust if needed)

REAL_OUT_DIR=~/PCP_genomics/KINETO/data/processed/phylogeny/histones/motifs/clustering/glam2_real_results
LOG_FILE=$REAL_OUT_DIR/../glam2_real_scores.log

# Start log file
echo "# File : Motif Score" > "$LOG_FILE"

# Loop through each GLAM2 result directory
for result_dir in "$REAL_OUT_DIR"/*; do
  if [ -d "$result_dir" ]; then
    base=$(basename "$result_dir")
    glam2_txt="$result_dir/glam2.txt"

    if [ -f "$glam2_txt" ]; then
      score=$(grep '^Score' "$glam2_txt" | awk '{print $2}')
      echo "$base : $score" >> "$LOG_FILE"
    else
      echo "Warning: glam2.txt not found in $result_dir"
    fi
  fi
done

echo "Motif scores saved to: $LOG_FILE"

