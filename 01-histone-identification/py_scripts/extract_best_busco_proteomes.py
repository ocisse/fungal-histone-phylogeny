import re
import argparse
def main():
    parser = argparse.ArgumentParser(description="Python script to get best scoring proteomes for each species")
    parser.add_argument("busco_batch_summary_file", help = "Batch summary file outputted by BUSCO")
    parser.add_argument("proteomes_for_hmmsearch", help = "Text file containing the names of the proteomes with the best BUSCO scores")
    args = parser.parse_args()
    parse_summary_file(args.busco_batch_summary_file, args.proteomes_for_hmmsearch)

def extract_species(proteome_file):
    """
    Extract species name from proteome filename, handling multiple formats:
    - GCA/GCF_NUMBER_Genus_species_strain_info_protein.faa
    - GCA/GCF_NUMBER_Genus species_strain_info_protein.faa
    """
    # Try first to match accession + full scientific name with space
    match = re.search(r'GC[AF]_\d+\.\d+_([A-Za-z]+\s+[A-Za-z]+)', proteome_file)
    if match:
        # For filenames with genus+species with a space between them
        return match.group(1)
    
    # Alternative match for underscore-separated genus and species
    match = re.search(r'GC[AF]_\d+\.\d+_([A-Za-z]+)_([A-Za-z]+)', proteome_file)
    if match:
        # For filenames with genus+species with an underscore between them
        return f"{match.group(1)} {match.group(2)}"
    
    # If all else fails, return the whole filename
    return proteome_file

def parse_summary_file(infile, outfile):
    with open(infile, "r") as fh_in:
        lines = fh_in.readlines()
    
    # Dictionary to store the highest scoring proteome for each species
    best_proteomes = {}
    
    # Debug counter
    line_count = 0
    
    for line in lines:
        line_count += 1
        line = line.strip()
        
        # Skip empty lines
        if not line:
            continue
            
        # Find position of eukaryota_odb10 to separate filename from metrics
        db_pos = line.find("eukaryota_odb10")
        if db_pos == -1:
            print(f"Warning: Line {line_count} doesn't contain expected database marker")
            continue
            
        # Extract filename and metrics
        proteome_file = line[:db_pos].strip()
        metrics_str = line[db_pos:].strip()
        metrics = metrics_str.split()
        
        # Extract species name
        species = extract_species(proteome_file)
        
        # Ensure we have enough metrics
        if len(metrics) >= 2:
            try:
                completeness = float(metrics[1])  # First value after db name
                # Store the best proteome for each species
                if species not in best_proteomes or completeness > best_proteomes[species][1]:
                    best_proteomes[species] = (proteome_file, completeness)
                    print(f"New best for {species}: {proteome_file} ({completeness}%)")
                else:
                    print(f"Skipping {proteome_file} for {species} ({completeness}% vs current best {best_proteomes[species][1]}%)")
            except ValueError as e:
                print(f"Error processing line {line_count}: {line}")
                print(f"Extracted metrics: {metrics}")
                print(f"Error: {e}")
    
    # Write results
    print(f"\nWriting {len(best_proteomes)} best proteomes to {outfile}")
    with open(outfile, "w") as fh_out:
        for species, (proteome_file, score) in sorted(best_proteomes.items()):
            fh_out.write(f"{proteome_file}\n")
            print(f"Best for {species}: {proteome_file} ({score}%)")

if __name__ == "__main__":
    main()
