#!/bin/bash

### Remove unnecessary header information from fasta file, leaving only ID

read -rp "Enter FASTA file path: " filepath
read -rp "Enter output file path: " outfile
grep -v "^>" "$filepath" | tr -d '\n' > "$outfile"
