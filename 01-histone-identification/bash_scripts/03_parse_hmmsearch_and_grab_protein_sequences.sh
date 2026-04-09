#!/bin/bash

# Bash script for parsing the tabular output produced by HMMsearch in order to extract the protein sequences
# Will use seqkit library/package for extracting the sequences from the proteomes

# Define required file paths 
phylogeny_dir="../../../data/processed/Phylogeny_Analysis"
histone_hmmsearch_output_dir="$phylogeny_dir/Histones/All_Fungal_Proteomes_Hmmsearch_Output"
proteomes_directory="../../../data/raw/unique_fungal_proteome_files"

# Create output directories for our sequences 

histone_hmmsearch_seqs_output_dir="$phylogeny_dir/Histones/All_Canonical_Histones_Sequences"
mkdir -p $histone_hmmsearch_seqs_output_dir

for file in $(ls $histone_hmmsearch_output_dir);do
  protein_ids=($(grep -v "#" $histone_hmmsearch_output_dir/$file | awk '{print $1}'))
  accession_id=$(basename $histone_hmmsearch_output_dir/$file | cut -d '_' -f 1,2)
  proteome_file=$(find "$proteomes_directory" -type f -name "${accession_id}*_protein.faa")

  for id in ${protein_ids[@]};do
   seqkit grep -rnp "$id" "$proteome_file" >> $histone_hmmsearch_seqs_output_dir/histone_seqs.fasta;
  done
done

