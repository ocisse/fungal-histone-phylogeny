from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse

def main():
    parser = argparse.ArgumentParser(description = "Testing for ideal cutoff for HmmSearch")
    parser.add_argument("sequence_file", help = "Sequence file")
    parser.add_argument("outputfile", help = "Truncated sequences")
    args = parser.parse_args()
    apply_truncation(args.sequence_file, args.outputfile)

def truncate_protein(seq, percent = 20):
    n_terminal = seq[:int(len(seq)/2)]
    truncate_length = int(0.1*len(n_terminal))
    truncated_n_terminal = n_terminal[truncate_length:]
    return truncated_n_terminal


def apply_truncation(seq_file, outfile):
    protein_records = [record for record in SeqIO.parse(seq_file, "fasta")]
    truncated_records = []
    for record in protein_records:
        truncated_seq = truncate_protein(record.seq)
        truncated_record = SeqRecord(truncated_seq, id=record.id)
        truncated_records.append(truncated_record)

    SeqIO.write(truncated_records, outfile, "fasta")


if __name__ == "__main__":
    main()
