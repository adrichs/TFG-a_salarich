#!/bin/bash

### Download protein sequences in fasta format from list of  NCBI IDs

# Prompt user for the output file name
read -p "Enter the name of the output file to save all sequences to: " output_file

# Prompt user for the input file name
read -p "Enter the name of the file containing the list of GI numbers: " input_file

# Check if input file exists
if [[ ! -f "$input_file" ]]; then
    echo "Error: File '$input_file' not found!"
    exit 1
fi

# Create/clear output file
> "$output_file"

# Read each GI number and download
while IFS= read -r gi_number; do
    # Skip empty lines
    [[ -z "$gi_number" ]] && continue
    
    echo "Downloading GI: $gi_number..."
    
    # Download from NCBI using efetch
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${gi_number}&rettype=fasta&retmode=text" | \
        awk -v id="$gi_number" '/^>/ {print ">GI_"id; next} {printf "%s", $0} END {print ""}' >> "$output_file"
    
    # Be polite to NCBI servers (they request max 3 requests/second)
    sleep 0.4
done < "$input_file"

echo "All sequences downloaded to $output_file"
