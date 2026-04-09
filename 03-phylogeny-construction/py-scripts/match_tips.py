import argparse
import csv
import re

def normalize_name(name):
    # Lowercase, remove punctuation, compress spaces
    name = name.lower()
    name = re.sub(r"[_\.]", " ", name)
    name = re.sub(r"\s+", " ", name)
    return name.strip()

def extract_species_from_tip(tip):
    # Remove the accession and suffix, leave species info
    parts = tip.split("_")
    if len(parts) < 3:
        return ""
    species_part = " ".join(parts[1:-1])  # Skip GCA and "protein"
    return normalize_name(species_part)

def extract_species_from_header(header):
    match = re.search(r"\[([^\[\]]+)\]", header)
    if match:
        return normalize_name(match.group(1))
    return ""

def main():
    parser = argparse.ArgumentParser(description="Find tips not matching any species in FASTA headers.")
    parser.add_argument("-t", "--tips", required=True, help="CSV file with 'Tip_Name' column")
    parser.add_argument("-f", "--fasta", required=True, help="FASTA file to search in headers")
    parser.add_argument("-o", "--output", default="unmatched_tips.csv", help="Output CSV for unmatched tips")
    args = parser.parse_args()

    # Load tip names and normalize
    with open(args.tips, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        tips = [row["Tip_Name"] for row in reader]
    tip_species = {tip: extract_species_from_tip(tip) for tip in tips}

    # Load FASTA headers and extract species names
    fasta_species_set = set()
    with open(args.fasta, "r") as f:
        for line in f:
            if line.startswith(">"):
                species = extract_species_from_header(line)
                if species:
                    fasta_species_set.add(species)

    # Match
    unmatched = []
    for tip, species in tip_species.items():
        if species not in fasta_species_set:
            unmatched.append(tip)

    # Output
    with open(args.output, "w", newline='') as out_csv:
        writer = csv.writer(out_csv)
        writer.writerow(["Unmatched_Tip_Name"])
        for tip in unmatched:
            writer.writerow([tip])

    print(f"Checked {len(tips)} tips. Found {len(unmatched)} unmatched. Saved to {args.output}")

if __name__ == "__main__":
    main()

