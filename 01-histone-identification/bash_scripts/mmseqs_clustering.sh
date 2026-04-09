#!/bin/bash 

# Must run python script first to create taxon id mapping for sequences 

# Bash script for clustering the sequences that were derived from Hmmsearch 
hmm_dir="AIP_2025/data/processed/HMM"
histone_fasta="$hmm_dir/Histones/Histone_seqs/histone_seqs.fasta"

tmp_dir_histones="$mmseqs2_histones_output/tmp_mmseq"
mmseqs2_histones_output="$hmm_dir/Histones/MMSEQS2_easyclust"

mmseqs easy-cluster $histone_fasta $mmseqs2_histones_output /tmp --cov-mode 0 -c 0.8 

# Next steps after mmseq clustering

#Create blast db for taxonomy classfication
tax_id_map="$hmm_dir/Histones/Histone_Taxonomy/taxon_id_to_seq_mapping.txt"
representative_histone_sequences="$hmm_dir/Histones/MMSEQS2_easyclust_rep_seq.fasta"

# Create query database from sequences
mmseqs createdb $representative_histone_sequences HistoneDB


# Create target database from NCNBI to search against 
mkdir ncbi-taxdump && cd ncbi-taxdump
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvzf taxdump.tar.gz

mmseqs createtaxdb HistoneDB /tmp --ncbi-tax-dump ncbi-tax-dump --tax-mapping-file $tax_id_map
mmseqs taxonomy HistoneDB HistoneDB taxonomyResult /tmp
mmseqs createtsv HistoneDB taxonomyResult taxonomyResult.tsv
mmseqs taxonomyreport HistoneDB taxonomyResult krona_report.html --report-mode 1 # Krona style report 
mmseqs taxonomyreport HistoneDB taxonomyResult report.html # Pavian
