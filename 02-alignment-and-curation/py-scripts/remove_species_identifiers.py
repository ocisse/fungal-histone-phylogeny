'''
Python script for removing identifier to fasta files
'''

import sys
from Bio import SeqIO

def main():
    fasta_file, fasta_minus_seq_identifiers_file = sys.argv[1:3]
    remove_seq_identifiers(fasta_file, fasta_minus_seq_identifiers_file)

def remove_seq_identifiers(fasta_file, outfile):
    with open(outfile, "w") as output_handle:
        records = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            record.id = ""
            if "_" in record.description:
                record.description = record.description.split("_", 1)[1]
            else:
                print(record.description)
            records.append(record)
        SeqIO.write(records, output_handle, "fasta")

if __name__ == "__main__":
    main()

