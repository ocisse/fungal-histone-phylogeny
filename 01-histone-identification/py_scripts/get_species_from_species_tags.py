import sys
from Bio import SeqIO, Phylo

def main():
    fasta_file, tree_file, tree_with_species_as_leafs = sys.argv[1:4]
    add_seq_identifiers(fasta_file, tree_file, tree_with_species_as_leafs)

def add_seq_identifiers(fasta_file, tree_file, outfile):
    species_tag_map = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        organism = record.description.split(']')[0].split('[')[1]
        species_tag = record.id.split('|')[0]
        species_tag_map[species_tag] = organism
     
    tree = Phylo.read(tree_file, "newick")
    for leaf in tree.get_terminals():
        species_tag = leaf.name.split('|')[0]
        leaf.name = f"{species_tag_map.get(species_tag)}" 
    Phylo.write(tree, outfile, "newick")

if __name__ == "__main__":
    main()


## Needed to manually remove double nested bracket species [[Candida] boidinii] for example                                                        
