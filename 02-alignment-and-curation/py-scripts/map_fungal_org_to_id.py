import ast
import sys 

def main():
    taxon_id_org_names_file, fasta_sequence_file, mapping_seq_to_id_outfile= sys.argv[1:4]
    mapping = map_taxon_id_to_rep_seqs(taxon_id_org_names_file, fasta_sequence_file)
    write_outfile(mapping_seq_to_id_outfile, mapping) 


def map_taxon_id_to_rep_seqs(taxon_id_file, rep_seq_file):
    '''
    Function to map the taxon id to each sequence within the fasta that was outputted
    by mmseqs2 clustering command
    '''
    with open(taxon_id_file, "r") as fh_in1:
        org_id_list = []
        organism_to_taxon_id = {}
        for line in fh_in1:
            tuple_data = ast.literal_eval(line.strip())
            org_id_list.append(tuple_data)
        for element in org_id_list:
            org_name = element[0]
            taxon_id = element[1]
            organism_to_taxon_id[org_name] = taxon_id

    with open(rep_seq_file, "r") as fh_in2:
        seq_id_to_taxon_id_mapping = {}
        for line in fh_in2:
            if line.startswith(">"):
                header = line.strip()
                sequence_accession = header.split()[0].split('>')[1]
                organism_name = header.split("[")[1].split("]")[0]
                if organism_name in organism_to_taxon_id:
                    seq_id_to_taxon_id_mapping[sequence_accession] = organism_to_taxon_id[organism_name]
    return seq_id_to_taxon_id_mapping


def write_outfile(mapping_seq_to_id_outfile, mapping):     
    with open(mapping_seq_to_id_outfile, "w") as fh_out:
        for accession_id, taxon_id in mapping.items():
            fh_out.write(f"{accession_id}\t{taxon_id}\n")


if __name__ == "__main__":
    main()

