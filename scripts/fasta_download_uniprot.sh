#!/bin/bash

### Download protein sequences in fasta format from list of  UniProt IDs

# Prompt user for the output file name and source file
read -p "Enter the name of the output file to save all sequences to: " output_file
read -p "Enter the name of the file containing the list of IDs: " input_file

# Check if input file exists
if [[ ! -f "$input_file" ]]; then
    echo "Error: File '$input_file' not found!"
    exit 1
fi

# Create/clear output file
> "$output_file"

# Read each accession number and download
while IFS= read -r accession; do
    # Skip empty lines
    [[ -z "$accession" ]] && continue
    
    echo "Downloading $accession..."
    curl -s "https://rest.uniprot.org/uniprotkb/${accession}.fasta" | \
        awk -v id="$accession" '/^>/ {print ">"id; next} {printf "%s", $0} END {print ""}' >> "$output_file"
    
    # Be polite to the server
    sleep 0.5
done < "$input_file"

echo "All sequences downloaded to $output_file"
