#!/bin/bash 

# Bash script for extracting only histone_h3 proteins from HMMsearch output

histones_dir="../../../data/processed/Phylogeny_Analysis/Histones"
hmmsearch_output_dir="$histones_dir/All_Fungal_Proteomes_Hmmsearch_Output"
mkdir -p $histones_dir/histone_h3_hmmsearch
outdir="$histones_dir/Histone_H3_Analysis/H3_HMMsearch_Output"


# Loop through the HMMsearch output files and grab all proteins that have scores >= 100 

for hmm_search_outfile in $(ls $hmmsearch_output_dir);do

    grep -v "^#" "$hmmsearch_output_dir/$hmm_search_outfile" | while read -r line; do
      score=$(echo $line | awk '{print $6}')

    if awk -v s="$score" 'BEGIN {exit !(s >= 100)}';then
      echo $line >> $outdir/$hmm_search_outfile
    fi;
    done < $hmmsearch_output_dir/$hmm_search_outfile
  done
