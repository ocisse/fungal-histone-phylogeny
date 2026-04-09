#python disorder_analysis.py ../Histone_Curated_Alignment_Unique_Ids_cenH3.fasta  --export-scores

python disorder_analysis.py Histone_Curated_Alignment_Unique_Ids_cenH3.fasta \
  --color-scheme publication \
  --format pdf \
  --dpi 600 \
  --font-family serif \
  --line-width 2.5 \
  --threshold 0.6 \
  --color-scheme pastel \
  --output-prefix cenh3_disorder

