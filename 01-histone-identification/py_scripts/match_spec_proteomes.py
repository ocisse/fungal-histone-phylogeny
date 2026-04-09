import re
import csv
import argparse
import difflib

def normalize_name(name):
    # Normalize species names for matching: lowercase, remove underscores, normalize spaces
    name = name.lower().replace("_", " ").replace("-", " ").strip()
    name = re.sub(r"\s+", " ", name)
    return name

def extract_species_from_header(header_line):
    match = re.search(r'>(Spec\d+)\|.*\[(.*?)\]', header_line)
    if match:
        spec_id = match.group(1)
        species_name = normalize_name(match.group(2))
        return spec_id, species_name
    return None, None

def build_species_dict(proteome_file):
    species_dict = {}
    with open(proteome_file) as f:
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) != 2:
                continue
            id_raw, size = parts
            parts2 = id_raw.split("_", 1)
            if len(parts2) < 2:
                continue
            species_part = normalize_name(parts2[1].replace("_protein", ""))
            species_dict[species_part] = size
    return species_dict

def main():
    parser = argparse.ArgumentParser(description="Map SpecX headers to proteome sizes.")
    parser.add_argument("-p", "--proteome", required=True, help="File 1: ID and proteome size")
    parser.add_argument("-f", "--fasta", required=True, help="File 2: FASTA headers with SpecX")
    parser.add_argument("-o", "--output", default="spec_proteome_map.csv", help="Output CSV")
    args = parser.parse_args()

    species_to_size = build_species_dict(args.proteome)
    species_keys = list(species_to_size.keys())

    output_rows = []
    with open(args.fasta) as f:
        for line in f:
            if line.startswith(">"):
                spec_id, fasta_species = extract_species_from_header(line)
                if spec_id and fasta_species:
                    # Try exact match first
                    if fasta_species in species_to_size:
                        size = species_to_size[fasta_species]
                        output_rows.append((spec_id, size))
                    else:
                        # Fuzzy match fallback
                        matches = difflib.get_close_matches(fasta_species, species_keys, n=1, cutoff=0.6)
                        if matches:
                            matched_species = matches[0]
                            size = species_to_size[matched_species]
                            output_rows.append((spec_id, size))
                            print(f"[Fuzzy match] '{fasta_species}' → '{matched_species}'")
                        else:
                            print(f"[Warning] No match for: '{fasta_species}'")

    # Write output CSV
    with open(args.output, "w", newline="") as out_file:
        writer = csv.writer(out_file)
        writer.writerow(["ID", "Proteome size"])
        writer.writerows(output_rows)

    print(f"Wrote {len(output_rows)} matched entries to {args.output}")

if __name__ == "__main__":
    main()

