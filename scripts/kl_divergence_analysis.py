#!/usr/bin/env python3
import numpy as np
import pandas as pd
from scipy.stats import entropy
import sys

### Calculate Kullback-Leibler divergence from the common positions of two sets of PSSMs using the PDB position identifiers

def kl_divergence(p, q, pseudocount=0.001):
    """Calculate KL divergence KL(P||Q) with pseudocounts to handle zeros"""
    p = np.array(p) + pseudocount
    q = np.array(q) + pseudocount
    p = p / p.sum()  # renormalize
    q = q / q.sum()
    return entropy(p, q)

def load_exposed_data(filepath):
    """Load exposed residue frequency data"""
    df = pd.read_csv(filepath, sep='\t')
    amino_acids = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    
    data = {}
    for _, row in df.iterrows():
        resnum = row['PDB_resnum']
        freqs = [row[aa] for aa in amino_acids]
        data[resnum] = freqs
    
    return data, amino_acids

def main():
    if len(sys.argv) != 3:
        print("Usage: python kl_divergence_analysis.py <thermostable_file> <generic_file>")
        sys.exit(1)
    
    thermo_file = sys.argv[1]
    generic_file = sys.argv[2]
    
    # Load data
    print("Loading thermostable data...")
    thermo_data, amino_acids = load_exposed_data(thermo_file)
    print("Loading generic data...")
    generic_data, _ = load_exposed_data(generic_file)
    
    # Find common positions
    common_positions = set(thermo_data.keys()) & set(generic_data.keys())
    print(f"\nFound {len(common_positions)} common positions")
    print(f"Thermostable positions: {len(thermo_data)}")
    print(f"Generic positions: {len(generic_data)}")
    
    # Calculate KL divergences
    results = []
    for pos in sorted(common_positions):
        thermo_freqs = thermo_data[pos]
        generic_freqs = generic_data[pos]
        
        # Calculate both directions
        kl_thermo_to_generic = kl_divergence(thermo_freqs, generic_freqs)
        kl_generic_to_thermo = kl_divergence(generic_freqs, thermo_freqs)
        
        # Average for symmetric measure
        kl_symmetric = (kl_thermo_to_generic + kl_generic_to_thermo) / 2
        
        results.append({
            'PDB_resnum': pos,
            'KL_thermo_to_generic': kl_thermo_to_generic,
            'KL_generic_to_thermo': kl_generic_to_thermo,
            'KL_symmetric': kl_symmetric,
            'thermo_entropy': entropy(np.array(thermo_freqs) + 0.001),  # Shannon entropy
            'generic_entropy': entropy(np.array(generic_freqs) + 0.001),
            'entropy_diff': entropy(np.array(thermo_freqs) + 0.001) - entropy(np.array(generic_freqs) + 0.001)
        })
    
    # Convert to DataFrame and sort by symmetric KL divergence
    df_results = pd.DataFrame(results)
    df_results = df_results.sort_values('KL_symmetric', ascending=False)
    
    # Add ranking
    df_results['rank'] = range(1, len(df_results) + 1)
    
    # Save results
    output_file = "kl_divergence_results.csv"
    df_results.to_csv(output_file, index=False)
    print(f"\nResults saved to {output_file}")
    
    # Print top 20 most discriminating positions
    print(f"\nTop 20 most discriminating positions:")
    print("Rank\tPosition\tKL_symmetric\tEntropy_diff")
    print("-" * 50)
    for _, row in df_results.head(20).iterrows():
        print(f"{row['rank']}\t{row['PDB_resnum']}\t{row['KL_symmetric']:.4f}\t\t{row['entropy_diff']:.4f}")

if __name__ == "__main__":
    main()
