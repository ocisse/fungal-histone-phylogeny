"""
Rules for calculating conservation scores from MSA
"""

rule calculate_conservation:
    input:
        msa = config["msa"]
    output:
        scores = f"{OUTDIR}/conservation/scores.txt",
        stats = f"{OUTDIR}/stats/conservation_summary.txt"
    params:
        method = config["conservation"]["method"],
        gap_threshold = config["conservation"]["gap_threshold"]
    threads: config["resources"]["conservation"]
    log:
        f"{OUTDIR}/logs/conservation.log"
      #conda:
  #  "../envs/conservation.yaml"
    script:
        "../../scripts/calculate_conservation.py"

rule plot_conservation:
    input:
        scores = rules.calculate_conservation.output.scores
    output:
        plot = f"{OUTDIR}/figures/conservation_profile.png"
    log:
        f"{OUTDIR}/logs/plot_conservation.log"
      #conda:
  #  "../envs/conservation.yaml"
    script:
        "../scripts/plot_conservation.py"
