'''
Script for downloading fungal proteome files from NCBI and changing names to human readable (from accession ID filenames to Species names)
'''

rule download_proteomes:
  shell:
    """mkdir -p ../../../../data/raw/all_fungal_proteomes
    cd ../../../../data/raw/all_fungal_proteomes && datasets download genome taxon fungi --include protein
    unzip ncbi_dataset.zip"""


rule rename_proteome_directories:
  input:
    json_file = "../../../../data/raw/all_fungal_proteomes/ncbi_dataset/assembly_data_report.json",
    base_directory = "../../../../data/raw/all_fungal_proteomes/ncbi_dataset/data"
  params:
      python_renaming_dirs_script = "change_fungal_proteome_dirname_from_accession_to_speciesname.py"
  shell:
    "python ../../python_scripts/{params.python_renaming_dir_script} {input.json_file} {input.base_directory} {input.json_file} {input.base_directory}"
