from Bio import AlignIO
import sys 

def main():
    infile = sys.argv[1]
    outfile = sys.argv[2]
    # Parse Stockholm file
    alignment = AlignIO.read(infile, "stockholm")

    # Extract sequence IDs
    seq_ids = [record.id.split("/")[0] for record in alignment]
    with open(outfile, "w") as fh_out:
        for id in seq_ids:
            fh_out.write(f"{id}\n")

if __name__ == "__main__":
    main()
