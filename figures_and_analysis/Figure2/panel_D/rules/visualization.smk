"""
Rules for generating ChimeraX visualization files
"""

rule identify_catd:
    input:
        pdb = rules.download_pdb.output.pdb,
        mapping = rules.map_conservation.output.mapping
    output:
        catd = f"{OUTDIR}/annotations/catd_region.txt"
    params:
        chain = config["pdb"]["chain"],
        approx_start = config["catd"]["start"],
        approx_end = config["catd"]["end"],
        motif = config["catd"]["search_motif"]
    log:
        f"{OUTDIR}/logs/catd.log"
      #conda:
  #  "../envs/conservation.yaml"
    script:
        "../../scripts/identify_catd.py"

rule conservation_statistics:
    input:
        mapping = rules.map_conservation.output.mapping,
        catd = rules.identify_catd.output.catd
    output:
        stats = f"{OUTDIR}/stats/conservation_stats.txt"
    log:
        f"{OUTDIR}/logs/stats.log"
      #conda:
  #  "../envs/conservation.yaml"
    script:
        "../../scripts/calculate_statistics.py"

rule generate_defattr:
    input:
        mapping = rules.map_conservation.output.mapping
    output:
        defattr = f"{OUTDIR}/chimerax/conservation.defattr"
    params:
        chain = config["pdb"]["chain"]
    log:
        f"{OUTDIR}/logs/defattr.log"
      #conda:
  #  "../envs/conservation.yaml"
    script:
        "../../scripts/generate_defattr.py"

rule generate_chimerax_script:
    input:
        defattr = rules.generate_defattr.output.defattr,
        catd = rules.identify_catd.output.catd,
        pdb = rules.download_pdb.output.pdb
    output:
        script = f"{OUTDIR}/chimerax/visualize.cxc"
    params:
        pdb_id = PDB_ID,
        chain = config["pdb"]["chain"],
        h4_chain = config["pdb"]["h4_chain"],
        dna_chains = config["pdb"]["dna_chains"],
        colors = config["visualization"]
    log:
        f"{OUTDIR}/logs/chimerax_script.log"
      #conda:
  #  "../envs/conservation.yaml"
    script:
        "../../scripts/generate_chimerax_script.py"

rule render_chimerax:
    input:
        script = rules.generate_chimerax_script.output.script,
        pdb = rules.download_pdb.output.pdb
    output:
        image = f"{OUTDIR}/figures/cnp1_conservation.png"
    params:
        outdir = f"{OUTDIR}/figures"
    log:
        f"{OUTDIR}/logs/render.log"
      #conda:
  #  "../envs/conservation.yaml"
    shell:
        """
        cd {params.outdir} && \
        chimerax --nogui --silent \
            --script "../../{input.script}" \
            2> ../../{log} || true
        """
