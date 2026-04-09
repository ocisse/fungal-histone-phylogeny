"""
Rules for aligning PDB structure to MSA
"""

rule download_pdb:
    output:
        pdb = f"{OUTDIR}/structures/{PDB_ID}.pdb"
    params:
        pdb_id = PDB_ID
    resources:
        downloads = 1
    log:
        f"{OUTDIR}/logs/download_{PDB_ID}.log"
    shell:
        """
        wget -O {output.pdb} \
            https://files.rcsb.org/download/{params.pdb_id}.pdb \
            2> {log}
        """

rule extract_pdb_sequence:
    input:
        pdb = f"{OUTDIR}/structures/{PDB_ID}.pdb"
    output:
        seq = f"{OUTDIR}/sequences/pdb_sequence.fasta"
    params:
        chain = config["pdb"]["chain"]
    log:
        f"{OUTDIR}/logs/extract_sequence.log"
      #conda:
  #  "../envs/conservation.yaml"
    script:
        "../../scripts/extract_pdb_sequence.py"

rule find_representative:
    input:
        msa = config["msa"]
    output:
        seq = f"{OUTDIR}/sequences/representative.fasta"
    params:
        keywords = config["keywords"]
    log:
        f"{OUTDIR}/logs/find_representative.log"
      #conda:
  #  "../envs/conservation.yaml"
    script:
        "../../scripts/find_representative.py"

rule align_pdb_to_msa:
    input:
        pdb_seq = f"{OUTDIR}/sequences/pdb_sequence.fasta",
        msa_seq = f"{OUTDIR}/sequences/representative.fasta"
    output:
        alignment = f"{OUTDIR}/alignments/pdb_to_msa.txt",
        stats = f"{OUTDIR}/stats/alignment_stats.txt"
    params:
        algorithm = config["alignment"]["algorithm"],
        min_identity = config["alignment"]["min_identity"]
    log:
        f"{OUTDIR}/logs/alignment.log"
      #conda:
  #  "../envs/conservation.yaml"
    script:
        "../../scripts/align_sequences.py"

rule map_conservation:
    input:
        alignment = f"{OUTDIR}/alignments/pdb_to_msa.txt",
        conservation = f"{OUTDIR}/conservation/scores.txt",
        pdb = f"{OUTDIR}/structures/{PDB_ID}.pdb"
    output:
        mapping = f"{OUTDIR}/mapping/pdb_to_conservation.txt",
        stats = f"{OUTDIR}/stats/mapping_stats.txt"
    params:
        chain = config["pdb"]["chain"]
    log:
        f"{OUTDIR}/logs/mapping.log"
      #conda:
  #  "../envs/conservation.yaml"
    script:
        "../../scripts/map_conservation.py"
