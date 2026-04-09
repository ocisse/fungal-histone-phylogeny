import argparse
from Bio import SeqIO, Phylo

def main():
    parser = argparse.ArgumentParser(description = "Script for adding species tags to Species Tree file")
    parser.add_argument("fasta_file", help = "Fasta file with species tags and organism names")
    parser.add_argument("tree_file", help = "Species Tree file")
    parser.add_argument("tree_plus_seq_identifiers_file", help = "Tree file with species tags")
    args = parser.parse_args()
    add_seq_identifiers(args.fasta_file, args.tree_file, args.tree_plus_seq_identifiers_file)

def add_seq_identifiers(fasta_file, tree_file, outfile):
    species_tag_map = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        organism = record.description.split(']')[0].split('[')[1]
        species_tag = record.id.split('|')[0]
        species_tag_map[organism] = species_tag
     
    tree = Phylo.read(tree_file, "newick")
    for leaf in tree.get_terminals():
        org_name = " ".join(leaf.name.split("_")[2:]).replace("_", " ").replace("protein", "").strip()
        leaf.name = f"{species_tag_map.get(org_name, 'NP')}" 
    Phylo.write(tree, outfile, "newick")

if __name__ == "__main__":
    main()


## Needed to manually remove double nested bracket species [[Candida] boidinii] for example                                                        
