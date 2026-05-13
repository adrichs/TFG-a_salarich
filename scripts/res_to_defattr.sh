#!/bin/bash

### Format .rsa to .defattr for ChimeraX compatibility

read -rp "Enter FreeSASA .rsa file path: " rsafile
read -rp "Enter output .defattr file path: " outfile

echo "attribute: relativeSASA
match mode: 1-to-1
recipient: residues" > "$outfile"

grep "^RES" "$rsafile" | while read -r _ resname chain resnum absall relall rest; do
    awk -v chain="$chain" -v resnum="$resnum" -v rel="$relall" \
        'BEGIN {printf "\t/%s:%s\t%.4f\n", chain, resnum, rel/100}'
done >> "$outfile"

echo "Done! Defattr file written to $outfile"
