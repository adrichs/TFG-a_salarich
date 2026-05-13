#!/bin/bash

### Remove identical IDs from two files to avoid dataset overlap

# Prompt for file names
read -p "Enter the name of the thermostable IDs file: " thermo_file
read -p "Enter the name of the generic IDs file: " generic_file
read -p "Enter the name for the output file (generic IDs without thermostable): " output_file

# Check if files exist
if [[ ! -f "$thermo_file" ]]; then
    echo "Error: File '$thermo_file' not found!"
    exit 1
fi

if [[ ! -f "$generic_file" ]]; then
    echo "Error: File '$generic_file' not found!"
    exit 1
fi

echo "Filtering IDs..."

# Create output file
> "$output_file"

# Read each ID from generic file
while IFS= read -r id; do
    # Skip empty lines
    [[ -z "$id" ]] && continue
    
    # Check if this ID exists in thermostable file
    if grep -qFx "$id" "$thermo_file"; then
        echo "  Removing: $id (found in thermostable set)"
    else
        echo "$id" >> "$output_file"
    fi
done < "$generic_file"

# Show statistics
thermo_count=$(grep -c . "$thermo_file")
generic_count=$(grep -c . "$generic_file")
output_count=$(grep -c . "$output_file")
removed_count=$((generic_count - output_count))

echo ""
echo "Summary:"
echo "  Thermostable IDs: $thermo_count"
echo "  Original generic IDs: $generic_count"
echo "  IDs removed: $removed_count"
echo "  Final generic IDs: $output_count"
echo ""
echo "Non-overlapping generic IDs saved to: $output_file"
