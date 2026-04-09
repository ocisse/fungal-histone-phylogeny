import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description = "Script to add unique species tags/identifiers to fasta files")
    parser.add_argument("fasta_file", help = "Original fasta file")
    parser.add_argument("fasta_plus_seq_identifiers_file", help = "Output file with unique species tags")
    args = parser.parse_args()
    add_seq_identifiers(args.fasta_file, args.fasta_plus_seq_identifiers_file)

def add_seq_identifiers(fasta_file, outfile):
    records = []
    unique_species = {}
    species_counter = 0

    for record in SeqIO.parse(fasta_file, "fasta"):
        org_name = record.description.split(']')[0].split('[')[1]

        if org_name not in unique_species:
            unique_species[org_name] = f"Spec{species_counter}"
            species_counter += 1
        
        record.id = f"{unique_species[org_name]}|"
        records.append(record)
    
    with open(outfile, "w") as out_handle:
        SeqIO.write(records, out_handle, "fasta")

if __name__ == "__main__":
    main()


## Needed to manually remove double nested bracket species [[Candida] boidinii] for example                                                        
