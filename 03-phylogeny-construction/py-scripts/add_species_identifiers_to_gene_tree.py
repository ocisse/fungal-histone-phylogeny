import argparse
import sys
from Bio import SeqIO, Phylo

def main():
    parser = argparse.ArgumentParser(description = "Script to add species tags to phylo trees, uses dictionary and maps to fasta files with organism names")
    parser.add_argument("fasta_file", help = "fasta file with species tags and organism names")
    parser.add_argument("tree_file", help = "gene tree file")
    parser.add_argument("tree_plus_seq_identifiers_file", help = "Gene tree file with organisms as leaf names for tax analysis")
    args = parser.parse_args()
    sys.setrecursionlimit(10000)
    add_seq_identifiers(args.fasta_file, args.tree_file, args.tree_plus_seq_identifiers_file)

def add_seq_identifiers(fasta_file, tree_file, outfile):
    species_tag_map = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        organism = record.description.split(']')[0].split('[')[1]
        species_tag = record.id.split('|')[0]
        species_tag_map[species_tag] = organism

    tree = Phylo.read(tree_file, "newick")

    for leaf in tree.get_terminals():
        species_tag = leaf.name.split("|")[0]
        leaf.name = f"{species_tag_map.get(species_tag)}"
    Phylo.write(tree, outfile, "newick")

if __name__ == "__main__":
    main()


## Needed to manually remove double nested bracket species [[Candida] boidinii] for example                                                        
