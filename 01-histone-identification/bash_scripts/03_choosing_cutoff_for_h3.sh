#!/bin/bash 

# Bash script for choosing the cutoff value for histone H3


# Download all swissprot proteins
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

histones_directory="../../../data/processed/Phylogeny_Analysis/Histones"
hmm_profile_histones="$histones_directory/hmm_profiles/PFAM_histone_HMM.hmm"
swissprot_database_fasta="$histones_directory/choosing_cutoff_for_histone_H3/uniprot_sprot.fasta"
swissprot_database="$histones_directory/choosing_cutoff_for_histone_H3/swissprot_database"
histone_H3_query="$histones_directory/choosing_cutoff_for_histone_H3/histone_H3_query_seq.fasta"
histone_H3_blastp_results="$histones_directory/choosing_cutoff_for_histone_H3/histone_H3_blast_results.out"

# Make blast db then run blastp using known Histone H3 sequence as query

makeblastdb -in $swissprot_database_fasta -dbtype prot -out $swissprot_database
blastp -query $histone_H3_query -db $swissprot_database -out histone_H3_blast_results.out -evalue 1e-5 -outfmt 6 -max_target_seqs 100

# Repeat for DNA polymerase, transcription factors, and non-histone histone-fold containing proteins
# This is to find precise cutoff score/threshold for choosing histone H3 from Hmmsearch output
# Must manually download the query sequences from Swissprot
dna_pol_query="$histones_directory/choosing_cutoff_for_histone_H3/dna_pol_query_seq.fasta"
dna_pol_blastp_results="$histones_directory/choosing_cutoff_for_histone_H3/dna_pol_blast_results.out"
blastp -query $dna_pol_query -db $swissprot_database -out $dna_pol_blastp_results -evalue 1e-5 -outfmt 6 -max_target_seqs 100


trans_factor_query="$histones_directory/choosing_cutoff_for_histone_H3/transcription_factor_query_seq.fasta"
tran_factor_blastp_results="$histones_directory/choosing_cutoff_for_histone_H3/trans_factor_blastp_results.out"
blastp -query $trans_factor_query -db $swissprot_database -out $trans_factor_blastp_results -evalue 1e-5 -outfmt 6 -max_target_seqs 100


histone_h4_query="$histones_directory/choosing_cutoff_for_histone_H3/histone_H4_query_seq.fasta"
histone_h4_blastp_results="$histones_directory/choosing_cutoff_for_histone_H3/histone_H4_blastp_results.out"
blastp -query $histone_h4_query -db $swissprot_database -out $histone_H4_blastp_results -evalue 1e-5 -outfmt 6 -max_target_seqs 100


histone_H2_query="$histones_directory/choosing_cutoff_for_histone_H3/histone_H2_query_seq.fasta"
histone_H2_blastp_results="$histones_directory/choosing_cutoff_for_histone_H3/histone_H4_blastp_results.out"
blastp -query $histone_H2_query -db $swissprot_database -out  -evalue 1e -5 -outfmt 6 -max_target_seqs 100

non_histone_hfd="$histones_directory/choosing_cutoff_for_histone_H3/non_histone_hfd_protein_seq.fasta"
non_histone_hfd_blastp_results="$histones_directory/choosing_cutoff_for_histone_H3/non_histone_hfd_blastp_results.out"
blastp -query $non_histone_hfd -db $swissprot_database -out $non_histone_hfd_blastp_results -evalue 1e-5 -outfmt 6 -max_target_seqs 100


