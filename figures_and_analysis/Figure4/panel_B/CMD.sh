python glam2_parser.py cdhit_cluster_fastas/cluster_10.fasta.glam2_out/glam2.txt --method first
python glam2_parser.py cdhit_cluster_fastas/cluster_64.fasta.glam2_out/glam2.txt --method first
python glam2_parser.py cdhit_cluster_fastas/cluster_56.fasta.glam2_out/glam2.txt --method first
python glam2_parser.py cdhit_cluster_fastas/cluster_306.fasta.glam2_out/glam2.txt --method first
python glam2_parser.py cdhit_cluster_fastas/cluster_131.fasta.glam2_out/glam2.txt --method first
python glam2_parser.py cdhit_cluster_fastas/cluster_38.fasta.glam2_out/glam2.txt --method first
python glam2_parser.py cdhit_cluster_fastas/cluster_98.fasta.glam2_out/glam2.txt --method first
python glam2_parser.py cdhit_cluster_fastas/cluster_91.fasta.glam2_out/glam2.txt --method first
python glam2_parser.py cdhit_cluster_fastas/cluster_14.fasta.glam2_out/glam2.txt --method first

head -2 cdhit_cluster_fastas/cluster_10.fasta.glam2_out/cluster_10_motifs_representatives.fasta > glam2_representatives.fasta
head -2 cdhit_cluster_fastas/cluster_64.fasta.glam2_out/cluster_64_motifs_representatives.fasta >> glam2_representatives.fasta
head -2 cdhit_cluster_fastas/cluster_56.fasta.glam2_out/cluster_56_motifs_representatives.fasta >> glam2_representatives.fasta 
head -2 cdhit_cluster_fastas/cluster_306.fasta.glam2_out/cluster_306_motifs_representatives.fasta >> glam2_representatives.fasta
head -2 cdhit_cluster_fastas/cluster_131.fasta.glam2_out/cluster_131_motifs_representatives.fasta >> glam2_representatives.fasta
head -2 cdhit_cluster_fastas/cluster_38.fasta.glam2_out/cluster_38_motifs_representatives.fasta >> glam2_representatives.fasta
head -2 cdhit_cluster_fastas/cluster_98.fasta.glam2_out/cluster_98_motifs_representatives.fasta >> glam2_representatives.fasta
head -2 cdhit_cluster_fastas/cluster_91.fasta.glam2_out/cluster_91_motifs_representatives.fasta >> glam2_representatives.fasta
head -2 cdhit_cluster_fastas/cluster_14.fasta.glam2_out/cluster_14_motifs_representatives.fasta >> glam2_representatives.fasta

python motif_ranker.py glam2_representatives.fasta --export-data 
python motif_ranker.py glam2_representatives.fasta --export-data --color-scheme publication --format pdf
