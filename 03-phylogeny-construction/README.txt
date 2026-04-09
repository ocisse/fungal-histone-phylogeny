#
## list of snakemake files

- 01_run_iqtree.smk # 
  rules:
  - "iqtree_HFD" -> compute ML phylogeny for HFD across all 6271 histone H3 like alignment
  - "iqtree_Nterm" -> compute ML phylogeny for N-term only 
  - "iqtree_Nterm_constrained" -> compute ML constrained phylogeny for N-term (pruned species tree as constraint)
  - "iqtree_Nterm_constrained_sub" -> same but constraint is sub tree sampled from species tree (50 tips)
  - "iqtree_Nterm_unconstrained" -> unconstrained ML phylogeny for N-term
  - "combine_nterm_cons_and_uncons_tree" -> topoplogy test; concatenate constrained and unconstrained trees
  - "Nterm_topology_testing" -> N-term tree topology testing
  - "iqtree_HFD_constrained" ->  compute ML constrained phylogeny for HFD (pruned species tree as constraint)
  - "iqtree_HFD_constrained_sub" -> same but subsampled species tree for HFD
  - "iqtree_HFD_unconstrained" ->  unconstrained ML phylogeny for HFD
  - "combine_hfd_cons_and_uncons_tree" ->  topoplogy test; concatenate constrained and unconstrained trees
  - "HFD_topology_testing" -> topoplogy testing
  - "combine_tree_full_with_wo_partition" -> combining ML tree without and with partition
  - "full_partition_topology_testing:"
