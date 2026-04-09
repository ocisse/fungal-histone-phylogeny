import argparse
import pandas as pd
import glob
import os

# Define the input directory containing your .txt files.
parser = argparse.ArgumentParser(description = "Python script to write orthologs to excel file")
parser.add_argument("input_directory", help = "directory where ortholog files are located")
parser.add_argument("ordered_proteins_file", help = "file with kinetochore proteins in correct order")
parser.add_argument("output_xlsx", help = "Outputted excel file ")
args = parser.parse_args()

# Create the file pattern for all .txt files.
file_pattern = os.path.join(args.input_directory, "*.txt")
file_list = glob.glob(file_pattern)


# Read the ordered protein names from the file into a list
with open(args.ordered_proteins_file, 'r') as f:
    protein_order = [line.strip() for line in f.readlines()]

# Dictionary to store data.
# Keys: protein identifiers (derived from the filename)
# Values: dictionaries mapping species name -> processed cell content.
protein_data = {}

# Set to collect all species names encountered.
all_species = set()

# Process each file.
for file_path in file_list:
    # Use the filename (without extension) as the protein identifier.
    protein_id = os.path.splitext(os.path.basename(file_path))[0]
    species_cells = {}
    
    # Add the protein_id to the order list (this keeps track of the order they appear)
    protein_order.append(protein_id)
    
    with open(file_path, 'r') as infile:
        lines = infile.readlines()
        
        # Ignore the first line (HOG ID)
        for line in lines[1:]:  # Start reading from the second line
            line = line.strip()
            if not line:
                continue  # Skip empty lines.
            
            # Split the line by comma to get all candidate proteins.
            candidates = [token.strip() for token in line.split(",")]
            
            # Use the first candidate to extract the species name.
            species = candidates[0].split("|", 1)[0].strip()
            
            # Process each candidate: extract the part after the pipe.
            protein_ids = []
            for candidate in candidates:
                parts = candidate.split("|", 1)
                if len(parts) == 2:
                    protein_ids.append(parts[1].strip())
                else:
                    # If there's no pipe, fall back to the full candidate.
                    protein_ids.append(candidate)
            
            # Join all protein IDs with a comma (or another separator) for the cell.
            cell_value = ", ".join(protein_ids)
            species_cells[species] = cell_value
            all_species.add(species)
    
    protein_data[protein_id] = species_cells

# Sort the species names to have a consistent column order.
species_list = sorted(all_species)

# Convert the dictionary to a pandas DataFrame, ensuring the protein order.
df = pd.DataFrame.from_dict(protein_data, orient="index", columns=species_list)

# Ensure the DataFrame follows the order of protein identifiers from the ordered file
df = df.loc[protein_order]

# Reset index so that "Protein" is the first column.
df.index.name = "Protein"
df.reset_index(inplace=True)

# Write the DataFrame to an Excel (.xlsx) file.
df.to_excel(args.output_xlsx, index=False, engine="openpyxl")


