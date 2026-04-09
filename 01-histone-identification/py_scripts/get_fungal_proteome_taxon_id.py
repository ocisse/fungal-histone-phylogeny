import argparse
import json


def main():
    parser = argparse.ArgumentParser(description = "Script to get the taxon ID of organisms from the json file downloaded along with proteomes using datasets")
    parser.add_argument("json_file", help = "json file with all proteome metadata downloaded via datasets")
    parser.add_argument("proteome_file", help = "proteome file")
    parser.add_argument("taxon_id_outfile", help = "outfile containing taxon IDs for each organism")
    args = parser.parse_args()
    accessions_dict = parse_json(args.json_file)
    taxon_ids = parse_proteome_file(args.proteome_file, accessions_dict)
    write_outfile(taxon_ids, args.taxon_id_outfile)


def parse_json(json_file):
    with open(json_file, "r") as fh_in:
        proteome_dict = [json.loads(line) for line in fh_in]
    #Build a dictionary mapping organism names to their taxonomic IDs
    accession_to_taxid = {}
    for entry in proteome_dict:
        accession = entry['accession']
        organism_name = entry['organism']['organismName']
        tax_id = entry['organism']['taxId']
        accession_to_taxid[accession] = organism_name, tax_id

    return accession_to_taxid


def parse_proteome_file(proteome_file, accession_to_taxid):
    with open(proteome_file, "r") as fh_in:
        taxon_ids = []
        for line in fh_in:
            accession_id = line.strip().split("_")[0:2]
            accession_id = "_".join(accession_id)
            if accession_id in accession_to_taxid:
                taxon_ids.append(accession_to_taxid[accession_id]) 
    return taxon_ids


def write_outfile(taxon_ids, outfile):
    with open(outfile, "w") as fh_out:
        for line in taxon_ids:
            fh_out.write(f"{line}\n")


if __name__ == "__main__":
    main()
