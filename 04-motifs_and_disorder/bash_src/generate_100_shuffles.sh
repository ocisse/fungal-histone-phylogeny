#!/bin/bash

# Exit if any command fails
set -e

# Define working directory
export DIR=~/PCP_genomics/KINETO/data/processed/phylogeny/histones/motifs/clustering

# Create output folder if it doesn't exist
mkdir -p "$DIR/shuffled_cdhit_cluster_fastas"

# Create or clear the seed log file
SEED_LOG="$DIR/shuffled_cdhit_cluster_fastas/shuffle_seeds.log"
echo "# File : Seed" > "$SEED_LOG"

# Loop over all FASTA files
for fasta in "$DIR/cdhit_cluster_fastas/"*.fasta; do
  base=$(basename "$fasta" .fasta)

  echo "Shuffling $base.fasta..."

  # Generate 100 shuffled versions
  for i in {1..100}; do
    seed=$RANDOM
    output_file="$DIR/shuffled_cdhit_cluster_fastas/${base}_shuffled_$i.fasta"

    # Shuffle the sequence
    fasta-shuffle-letters -kmer 1 -protein -seed "$seed" "$fasta" "$output_file"

    # Log the filename and seed used
    echo "${base}_shuffled_$i.fasta : seed=$seed" >> "$SEED_LOG"
  done
done

echo "Done generating 100 shuffled FASTAs per input file."
echo "Seeds logged in: $SEED_LOG"

