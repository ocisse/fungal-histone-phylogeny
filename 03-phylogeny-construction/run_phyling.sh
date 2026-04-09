#!/bin/bash 

#SBATCH --job-name=phyling
#SBATCH --partition=norm
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=10g
#SBATCH --time=5-00:00


# Bash script for running PHYling to create Species Tree

# First clone the PHYling github repo and create the conda env from the repo's env.yaml file
# git clone https://github.com/stajichlab/PHYling.git
# pip install .

# Now ready for analysis 
eval "$(mamba shell hook --shell bash)"
mamba activate phyling


unique_fungal_proteome_files="../../../data/raw/unique_fungal_proteome_files"

phyling_align_output="../../../data/processed/PHYling/PHYling_align_output"
busco_eukaryota_db_dir="../../../data/processed/BUSCO_data/busco_downloads/lineages/eukaryota_odb10/hmms"
phyling align -I $unique_fungal_proteome_files -o $phyling_align_output -m $busco_eukaryota_db_dir

phyling_filtering_output="../../../data/processed/PHYling/PHYling_filtering_output"
mkdir -p $phyling_filtering_output
phyling filter -I $output_dir -n 20 -t 12 -o $phyling_filtering_output 

# Create tree with IQ-Tree
phyling_tree_output="../../../data/processed/PHYling/PHYling_iqtree_output"
mkdir -p $phyling_tree_output
# phyling tree -I $phyling_filtering_output -M iqtree -t 16 -f -c -v -p -o $phyling_tree_output
iqtree2 -s $phyling_tree_output/concat_alignments.mfa --prefix \
  $phyling_tree_output/iqtree/concat_alignments.mfa -T AUTO --threads-max 16 \
  --seqtype AA -p $phyling_tree_output/concat_alignments.partition #--redo-tree
