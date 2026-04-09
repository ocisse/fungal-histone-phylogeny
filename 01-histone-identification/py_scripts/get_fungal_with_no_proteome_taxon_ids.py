import json
import sys

def main():
    json_file, taxon_id_outfile = sys.argv[1:3]
    taxon_ids = parse_json(json_file)
    write_outfile(taxon_ids, taxon_id_outfile)


def parse_json(json_file):
    with open(json_file, "r") as fh_in:
        proteome_dict = [json.loads(line) for line in fh_in]
    #Build a dictionary mapping organism names to their taxonomic IDs
    accession_to_taxid = {}
    tax_ids = []
    for entry in proteome_dict:
        accession = entry['accession']
        tax_id = entry['organism']['tax_id']
        accession_to_taxid[accession] = tax_id
        tax_ids.append(tax_id)
    return tax_ids


def write_outfile(tax_ids, outfile):
    with open(outfile, "w") as fh_out:
        for line in tax_ids:
            fh_out.write(f"{line}\n")


if __name__ == "__main__":
    main()
