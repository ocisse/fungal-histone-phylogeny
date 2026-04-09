import argparse
from pprint import pprint 
def main():
    parser = argparse.ArgumentParser(description = "Script to create dictionary of all orthogroups")
    parser.add_argument("infile", help = "TSV file of orthogroups outputted by Orthofinder")
    parser.add_argument("outfile", help = "Orthogroup dictionary outfile")
    args = parser.parse_args()

    with open(args.infile, "r") as handle:
        ortho_dict = {}

        for line in handle:
            line = [x for x in line.rstrip().split("\t") if x.strip()]
            group = line[0]
            genes = line[3:]
            ortho_dict[group] = genes
    
    with open(args.outfile, "w") as handle:
        pprint(ortho_dict, stream = handle)

if __name__ == "__main__":
    main()
