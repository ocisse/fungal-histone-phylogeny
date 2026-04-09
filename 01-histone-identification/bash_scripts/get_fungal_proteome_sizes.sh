#!/bin/bash

# Bash script for getting the proteome sizes for all fungal species downloaded from NCBI

proteomes_dir="AIP_2025/data/raw/all_fungal_proteomes/ncbi_dataset/data"
output_dir="/data/littlestoneed/AIP_2025/data/processed/BUSCO_data"

#!/bin/bash

proteomes_dir="AIP_2025/data/raw/all_fungal_proteomes/ncbi_dataset/data"
output_dir="AIP_2025/data/processed/BUSCO_data"

for species_dir in "$proteomes_dir"/*/; do
  if [[ -d "$species_dir" ]]; then
    species_name=$(basename "$species_dir")
    proteome_file="$species_dir/protein.faa"
    if [[ -f "$proteome_file" ]]; then
      count=$(grep -c ">" "$proteome_file")
      printf "%s\t%s\n" "$species_name" "$count" >> "$output_dir/proteome_sizes.txt"
    fi
  fi
done

