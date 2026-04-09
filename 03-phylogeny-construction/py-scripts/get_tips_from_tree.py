import argparse
import csv
from Bio import Phylo

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Extract all tip names from a phylogenetic tree and save to CSV.")
    parser.add_argument("-i", "--input", required=True, help="Input tree file (e.g., tree.nwk)")
    parser.add_argument("-f", "--format", default="newick", help="Tree file format (default: newick)")
    parser.add_argument("-o", "--output", default="tree_tips.csv", help="Output CSV file (default: tree_tips.csv)")
    
    args = parser.parse_args()

    # Load the tree
    tree = Phylo.read(args.input, args.format)

    # Extract tip names
    tip_names = [term.name for term in tree.get_terminals()]

    # Write to CSV
    with open(args.output, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Tip_Name"])  # Header
        for name in tip_names:
            writer.writerow([name])

    print(f"Successfully wrote {len(tip_names)} tips to {args.output}")

if __name__ == "__main__":
    main()

