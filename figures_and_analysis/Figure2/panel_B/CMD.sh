CIAlign --infile Histone_Curated_Alignment_Unique_Ids_HFD.fasta \
  --outfile HFD --plot_consensus_similarity --plot_format svg --plot_similarity_palette viridis

python alignment_identity.py Histone_Curated_Alignment_Unique_Ids_NTERM_revised.fasta \
  --include-gaps \
  --output Histone_Curated_Alignment_Unique_Ids_NTERM_revised.identity_stats.txt \
  --plot Histone_Curated_Alignment_Unique_Ids_NTERM_revised.identity_distribution.png
