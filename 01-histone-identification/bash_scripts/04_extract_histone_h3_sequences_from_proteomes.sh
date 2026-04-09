#!/bin/bash

# Bash script for extracting histone h3 sequences from proteomes

histone_h3_hmmsearch_output="../../../data/processed/Phylogeny_Analysis/Histones/Histone_H3_Analysis/H3_HMMsearch_Output"
proteome_directory="../../../data/raw/unique_fungal_proteome_files"
histone_h3_sequence_output="../../../data/processed/Phylogeny_Analysis/Histones/Histone_H3_Analysis/H3_Sequences/histone_h3_sequences.fasta"


for file in $(ls $histone_h3_hmmsearch_output);do
  protein_ids=($(grep -v "#" $histone_h3_hmmsearch_output/$file | awk '{print $1}'))
  accession_id=$(basename $histone_h3_hmmsearch_output/$file | cut -d '_' -f 1,2)
  proteome_file=$(find "$proteome_directory" -type f -name "${accession_id}*_protein.faa")

  tmp_id_file=$(mktemp)
  printf "%s\n" "${protein_ids[@]}" > "$tmp_id_file"
  seqtk subseq "$proteome_file" "$tmp_id_file" >> $histone_h3_sequence_output
  rm "$tmp_id_file"
done
