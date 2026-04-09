#!/bin/bash

## Bash script for grabbing sequences from blastp results 

phylogeny_analysis="../../../data/processed/Phylogeny_Analysis"
h3_cutoff_dir="$phylogeny_analysis/Histones/choosing_cutoff_for_histone_H3"
blastp_results=(histone_H2 histone_H3 histone_H4 dna_pol non_histone_hfd)

for protein in ${blastp_results[@]};do
  while read -r line; do 
    protein_ids=($(awk '{print $2}' | cut -d '|' -f 2))
    tmp_id_file=$(mktemp)

    printf "%s\n" "${protein_ids[@]}" > "$tmp_id_file"
    for id in $tmp_id_file; do
    seqkit grep -rnp "$id" >> $h3_cutoff_dir/${protein}_blastp_sequences.fasta;
  done
  rm $tmp_id_file
  done<$h3_cutoff_dir/${protein}_blastp_results.out;
done


hmmsearch 
