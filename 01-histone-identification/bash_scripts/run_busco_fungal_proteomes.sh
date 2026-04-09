# Run BUSCO on proteomes to assess completedness
busco_outdir="/data/littlestoneed/AIP_2025/data/processed/BUSCO_data/Fungal_genomes_Busco_scores"
busco -i $data_dir/proteomes_for_busco \
  -m proteins \
  -l eukaryota_odb10 --cpu $SLURM_CPUS_PER_TASK \
  --out $busco_outdir \
  -r \
  --cpu $SLURM_CPUS_PER_TASK
