# Code and Analysis for: [Your Paper Title Here]

[![Paper DOI](https://img.shields.io/badge/Paper%20DOI-10.XXXX/journal.XXXXXX-blue)](https://doi.org/10.XXXX/journal.XXXXXX)
[![Zenodo DOI](https://img.shields.io/badge/Zenodo-10.5281/zenodo.XXXXXXX-blue)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains all the code, scripts, and workflows used to perform the large-scale phylogenetic analysis of fungal histones described in [**Your Last Name et al., Journal, Year**].

The primary goal of this project was to analyze over 6,000 fungal proteomes to identify histone candidates, curate them, and build a robust large-scale phylogeny to investigate their evolutionary history, with a focus on CENP-A homologs.

---

## Table of Contents
- [Setup and Dependencies](#setup-and-dependencies)
- [Data Availability](#data-availability)
- [Analysis Workflow](#analysis-workflow)
  - [Step 1: Histone Identification and Proteome Curation](#step-1-histone-identification-and-proteome-curation)
  - [Step 2: Alignment and Curation](#step-2-alignment-and-curation)
  - [Step 3: Phylogeny Construction](#step-3-phylogeny-construction)
  - [Step 4: Motif Discovery and Disorder Prediction](#step-4-motif-discovery-and-disorder-prediction)
- [Figure Generation](#figure-generation)
- [Repository Structure](#repository-structure)
- [Citation](#citation)
- [License](#license)

---

## Setup and Dependencies

All computational analyses were performed within a Conda environment. To ensure reproducibility, all required software and dependencies are specified in the `environment.yml` file.

1.  **Create the Conda environment:**
    ```bash
    conda env create -f environment.yml
    ```

2.  **Activate the environment before running any scripts:**
    ```bash
    conda activate fungal-histone-phylogeny
    ```
    *(Note: The environment name is defined inside the `environment.yml` file.)*

---

## Data Availability

Due to their size, all large data files are not stored in this GitHub repository. This includes:
- Raw fungal proteomes
- HMM profiles
- Large-scale sequence alignments
- Intermediate and final phylogenetic trees (`.treefile`, `.tre`, etc.)
- Large data matrices and other outputs

All these datasets are permanently archived on **Zenodo**. The directory structure in the Zenodo archive mirrors the structure of this GitHub repository.

**Zenodo Archive DOI:** [10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX)

To run the workflows, you must first download and place the data into the appropriate directories within this repository structure.

---

## Analysis Workflow

The analysis is structured as a series of sequential stages, with each stage located in a numbered directory. The primary data processing pipelines are orchestrated using a combination of Snakemake and bash scripts. **These workflows should be run in order.**

### Step 1: Histone Identification and Proteome Curation (`01-histone-identification/`)

This stage handles the initial download of >6000 fungal proteomes, quality control using BUSCO, and identification of histone candidates using HMM searches.

- **Key scripts:**
  - `01_Snakefile_download_fungal_proteomes_from_ncbi`: Downloads proteomes.
  - `02_Snakefile_busco_on_all_proteomes`: Runs BUSCO for quality control.
  - `bash_scripts/02_run_hmmsearch.sh`: Runs HMMsearch to find histone candidates.
- **To execute:** Follow the instructions within the `01-histone-identification/` directory. Each Snakefile can be run independently (e.g., `snakemake -s 01_Snakefile...`).

### Step 2: Alignment and Curation (`02-alignment-and-curation/`)

This stage takes the candidate sequences from Step 1, performs multiple sequence alignment (MSA) using MAFFT/MUSCLE, and cleans the alignment using tools like CIAlign.

- **Key scripts:**
  - `Snakefile_run_muscle`: Orchestrates the MSA.
  - `py-scripts/clean_alignment.py`: Custom scripts for refining the alignment.
- **To execute:** The main workflow is `Snakefile_run_muscle`.

### Step 3: Phylogeny Construction (`03-phylogeny-construction/`)

This stage uses the final, curated alignment to build the large-scale phylogenies with IQ-TREE and BEAST. It also includes workflows for topology testing and site-rate evolution analysis.

- **Key workflows:**
  - `workflow/01_Snakefile_iqtree3_histone...`: Runs IQ-TREE for ML phylogeny.
  - `run_beast.sh`: Script for running BEAST for Bayesian phylogeny.
  - `workflow/02_Snakefile_Topology_testing.sn`: Performs topology tests.

### Step 4: Motif Discovery and Disorder Prediction (`04-motifs_and_disorder/`)

This stage analyzes the curated sequences for conserved motifs (using GLAM2) and predicts intrinsically disordered regions (using IUPRED3).

- **Key workflow:** `Snakefile_motif_discovery`
- **To execute:** Run `snakemake -s Snakefile_motif_discovery`.

---

## Figure Generation

All scripts used to generate the final figures for the manuscript are located in the `figures_and_analysis/` directory. Each figure has its own subdirectory containing the Quarto notebook (`.qmd`), code, and necessary metadata to reproduce the plot.

To regenerate a specific figure, navigate to its directory and render the Quarto document. For example, for Figure 1:
```bash
cd figures_and_analysis/Figure1
quarto render fig1.qmd
