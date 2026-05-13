#!/usr/bin/env python3

### Extract exposed or buried amino acids from the PSSMs generated with Goalign, using a list of residues obtained with ChimeraX

# Three-letter to one-letter amino acid conversion
THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def read_pdb_sequence(pdb_file, chain='A'):
    """
    Read actual residue numbers and sequence from PDB ATOM records.
    Returns a list of (resnum, one_letter_code) tuples in order,
    with each residue appearing only once.
    """
    seen = set()
    residues = []
    with open(pdb_file) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            chain_id = line[21]
            if chain_id != chain:
                continue
            resname = line[17:20].strip()
            resnum = int(line[22:26].strip())
            if resnum in seen:
                continue
            seen.add(resnum)
            one_letter = THREE_TO_ONE.get(resname)
            if one_letter is None:
                continue
            residues.append((resnum, one_letter))
    return residues

def read_fasta_sequence(fasta_file, identifier='P07788'):
    """
    Read sequence of a specific entry from a FASTA file.
    Returns the full sequence string including gap characters.
    """
    seq = ''
    found = False
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>' + identifier):
                found = True
                continue
            if found and line.startswith('>'):
                break
            if found:
                seq += line.upper()
    return seq

def remove_gaps(seq):
    """Remove gap characters from a sequence string."""
    return ''.join(aa for aa in seq if aa not in ('-', '.'))

def read_exposed_residues(exposed_file):
    """Read PDB residue numbers from ChimeraX output file."""
    exposed = set()
    with open(exposed_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            resnum = int(parts[2].split(':')[1])
            exposed.add(resnum)
    return exposed

def read_pssm(pssm_file):
    """Read full PSSM file into a dictionary keyed by row number."""
    amino_acids = []
    pssm = {}
    with open(pssm_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if not parts[0].isdigit():
                amino_acids = [p for p in parts if p.strip() in list('ARNDCQEGHILKMFPSTWYV')]
                continue
            pos = int(parts[0])
            freqs = {aa: float(v) for aa, v in zip(amino_acids, parts[1:])}
            pssm[pos] = freqs
    return pssm, amino_acids

def build_canonical_to_pdb_mapping(orig_nogap, pdb_residues):
    """
    Map canonical sequence positions (1-based) to real PDB residue numbers.
    Walks both sequences simultaneously, matching amino acids.
    Positions in orig_nogap that are absent from the PDB (disordered) get no mapping.
    """
    canonical_to_pdb = {}
    pdb_idx = 0
    for orig_idx, aa in enumerate(orig_nogap):
        canonical_pos = orig_idx + 1
        if pdb_idx < len(pdb_residues) and aa == pdb_residues[pdb_idx][1]:
            canonical_to_pdb[canonical_pos] = pdb_residues[pdb_idx][0]
            pdb_idx += 1
    return canonical_to_pdb

def build_trim_pos_to_pssm_row(trim_seq):
    """
    Map each non-gap position in the trimmed alignment sequence (0-indexed)
    to its actual PSSM row number (1-indexed alignment column).
    This is necessary because P07788 may have gaps in some alignment columns,
    meaning PSSM row numbers do not equal gap-free position indices.
    """
    trim_pos_to_pssm = {}
    trim_pos = 0
    for col_idx, aa in enumerate(trim_seq):
        if aa not in ('-', '.'):
            trim_pos_to_pssm[trim_pos] = col_idx + 1
            trim_pos += 1
    return trim_pos_to_pssm

def build_pdb_to_pssm_mapping(orig_nogap, trim_nogap, trim_seq, canonical_to_pdb):
    """
    Map real PDB residue numbers to PSSM row numbers.
    Uses the actual alignment column index as the PSSM row number,
    correctly handling cases where P07788 has gaps in the trimmed alignment.
    """
    trim_pos_to_pssm = build_trim_pos_to_pssm_row(trim_seq)

    pdb_to_pssm = {}
    trim_idx = 0
    for orig_idx, aa in enumerate(orig_nogap):
        canonical_pos = orig_idx + 1
        if trim_idx < len(trim_nogap) and aa == trim_nogap[trim_idx]:
            pssm_row = trim_pos_to_pssm[trim_idx]  # actual column number
            trim_idx += 1
            if canonical_pos in canonical_to_pdb:
                real_pdb_resnum = canonical_to_pdb[canonical_pos]
                pdb_to_pssm[real_pdb_resnum] = pssm_row
    return pdb_to_pssm

def main():
    print("=== Exposed/Buried PSSM Extractor ===\n")
    orig_msa    = input("Enter original (untrimmed) MSA file path: ").strip()
    trim_msa    = input("Enter trimmed MSA file path: ").strip()
    pssm_file   = input("Enter full PSSM file path (from goalign): ").strip()
    exposed_file = input("Enter exposed residues file path (ChimeraX output): ").strip()
    pdb_file    = input("Enter PDB file path (e.g. 1GSK.pdb): ").strip()
    outfile     = input("Enter output file path: ").strip()

    print("\nReading PDB structure...")
    pdb_residues = read_pdb_sequence(pdb_file, chain='A')
    print(f"  {len(pdb_residues)} residues with coordinates found in PDB chain A")

    print("Reading exposed residues...")
    exposed = read_exposed_residues(exposed_file)
    print(f"  {len(exposed)} exposed/buried positions found")

    print("Reading original MSA...")
    orig_seq = read_fasta_sequence(orig_msa, identifier='P07788')
    orig_nogap = remove_gaps(orig_seq)
    print(f"  P07788 canonical sequence length: {len(orig_nogap)} residues")

    print("Reading trimmed MSA...")
    trim_seq = read_fasta_sequence(trim_msa, identifier='P07788')
    trim_nogap = remove_gaps(trim_seq)
    print(f"  P07788 trimmed sequence length: {len(trim_nogap)} residues")
    print(f"  Positions removed by trimAl: {len(orig_nogap) - len(trim_nogap)}")

    print("Reading PSSM...")
    pssm, amino_acids = read_pssm(pssm_file)
    print(f"  {len(pssm)} PSSM rows read")

    print("Building canonical sequence -> PDB residue number mapping...")
    canonical_to_pdb = build_canonical_to_pdb_mapping(orig_nogap, pdb_residues)
    print(f"  {len(canonical_to_pdb)} canonical positions mapped to real PDB numbers")
    print(f"  {len(orig_nogap) - len(canonical_to_pdb)} canonical positions absent from PDB (disordered)")

    print("Building PDB residue number -> PSSM row mapping...")
    pdb_to_pssm = build_pdb_to_pssm_mapping(orig_nogap, trim_nogap, trim_seq, canonical_to_pdb)
    print(f"  {len(pdb_to_pssm)} PDB positions mapped to PSSM rows")

    print(f"Writing output to {outfile}...")
    skipped = 0
    written = 0
    with open(outfile, 'w') as out:
        out.write('PDB_resnum\tPSSM_row\t' + '\t'.join(amino_acids) + '\n')
        for resnum in sorted(exposed):
            if resnum not in pdb_to_pssm:
                skipped += 1
                continue
            pssm_row = pdb_to_pssm[resnum]
            if pssm_row not in pssm:
                skipped += 1
                continue
            freqs = pssm[pssm_row]
            out.write(
                str(resnum) + '\t' +
                str(pssm_row) + '\t' +
                '\t'.join(f"{freqs[aa]:.3f}" for aa in amino_acids) + '\n'
            )
            written += 1

    print(f"\nDone! {written} positions written to {outfile}")
    print(f"Skipped {skipped} positions (trimmed by trimAl or absent from PDB)")

if __name__ == "__main__":
    main()
