#!/bin/bash

# Run hmmsearch with PFAM HMM Profiles to get our proteins of interest from our proteomes

# Setup variables for directory/file paths

proteomes_directory="../../../data/raw/unique_fungal_proteome_files"
phylogeny_dir="../../../data/processed/Phylogeny_Analysis"
mkdir -p $phylogeny_dir/Fungal_Proteomes_for_Hmmsearch


pfam_hmm_histone_profile="$phylogeny_dir/Histones/hmm_profiles/PFAM_histone_HMM.hmm"

histone_output_dir="$phylogeny_dir/Histones/All_Fungal_Proteomes_Hmmsearch_Output"

# First run hmmsearch to extract all histones
for proteome in $(ls $proteomes_directory); do
  output_file="$histone_output_dir/$(basename $proteome)"
  hmmsearch -T 23 --tblout "$output_file" "$pfam_hmm_histone_profile" "$proteomes_directory/$proteome";
done


