"""
Script for changing the directory names of the proteomes downloaded from NCBI from the accession ID to the Species Name
"""

import argparse
import json
import os

def main():
    parser = argparse.ArgumentParser(description="Change directory names from accession ID to Species Names")
    parser.add_argument("json_file", help = "Path to the ncbi json file containing metadata for all proteomes")
    parser.add_argument("base_dir", help = "base directory where proteome directories are located")
    args = parser.parse_args()
    #json_file, base_dir = sys.argv[1:3]
    accessions_dict = parse_json(args.json_file)
    rename_dirs(args.base_dir, accessions_dict)    

def parse_json(json_file):
    with open(json_file, "r") as fh_in:
        genome_dict = [json.loads(line) for line in fh_in]

    organism_accessions = {}
    for entry in genome_dict:
        accession = entry['accession']
        organism_name = entry['organism']['organismName']
        organism_accessions[accession] = organism_name
    
    return organism_accessions

def rename_dirs(base_dir, organism_accessions):
    for folder in os.listdir(base_dir):
        safe_organism_name = organism_accessions.get(folder, "Unknown").replace("/", "_").replace(" ", "_")
        new_name = f"{folder}_{safe_organism_name}"
        old_path = os.path.join(base_dir, folder)
        new_path = os.path.join(base_dir, new_name)
        os.rename(old_path, new_path)

if __name__ == "__main__":
    main()
