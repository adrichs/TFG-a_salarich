#!/usr/bin/env python3
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact
from collections import defaultdict
import sys

### Calculate amino acid group enrichment, using IMGT classification

# Define amino acid groups by properties (official IMGT classifications)
AA_GROUPS = {
    # IMGT Hydropathy classes
    'Hydrophobic': ['A', 'C', 'I', 'L', 'M', 'F', 'W', 'V'],
    'Neutral': ['G', 'H', 'P', 'S', 'T', 'Y'],
    'Hydrophilic': ['R', 'N', 'D', 'Q', 'E', 'K'],
    
    # IMGT Volume classes
    'Very_small': ['A', 'G', 'S'],
    'Small': ['N', 'D', 'C', 'P', 'T'],
    'Medium': ['Q', 'E', 'H', 'V'],
    'Large': ['R', 'I', 'L', 'K', 'M'],
    'Very_large': ['F', 'W', 'Y'],
    
    # IMGT Chemical classes
    'Aliphatic': ['A', 'G', 'I', 'L', 'P', 'V'],
    'Aromatic': ['F', 'W', 'Y'],
    'Sulfur': ['C', 'M'],
    'Hydroxyl': ['S', 'T'],
    'Basic': ['R', 'H', 'K'],
    'Acidic': ['D', 'E'],
    'Amide': ['N', 'Q'],
    
    # IMGT Charge classes
    'Positive_charged': ['R', 'H', 'K'],
    'Negative_charged': ['D', 'E'],
    'Uncharged': ['A', 'N', 'C', 'Q', 'G', 'I', 'L', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'],
    
    # IMGT Polarity classes
    'Polar': ['R', 'N', 'D', 'Q', 'E', 'H', 'K', 'S', 'T', 'Y'],
    'Nonpolar': ['A', 'C', 'G', 'I', 'L', 'M', 'F', 'P', 'W', 'V'],

    # Combined charge class
    'Charged': ['R', 'H', 'K', 'D', 'E']
}

def load_exposed_data(filepath):
    """Load exposed residue frequency data"""
    df = pd.read_csv(filepath, sep='\t')
    amino_acids = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    
    data = {}
    for _, row in df.iterrows():
        resnum = row['PDB_resnum']
        freqs = [row[aa] for aa in amino_acids]
        data[resnum] = dict(zip(amino_acids, freqs))
    
    return data, amino_acids

def calculate_group_enrichment(thermo_data, generic_data, common_positions, group_name, group_aas):
    """Calculate enrichment for a specific amino acid group"""
    
    thermo_group_total = 0
    thermo_other_total = 0
    generic_group_total = 0
    generic_other_total = 0
    
    for pos in common_positions:
        thermo_freqs = thermo_data[pos]
        generic_freqs = generic_data[pos]
        
        # Normalize frequencies to handle gaps
        thermo_sum = sum(thermo_freqs.values())
        generic_sum = sum(generic_freqs.values())
        
        if thermo_sum > 0:
            thermo_freqs_norm = {aa: freq/thermo_sum for aa, freq in thermo_freqs.items()}
        else:
            thermo_freqs_norm = thermo_freqs
            
        if generic_sum > 0:
            generic_freqs_norm = {aa: freq/generic_sum for aa, freq in generic_freqs.items()}
        else:
            generic_freqs_norm = generic_freqs
        
        # Sum frequencies for this group vs others (using normalized frequencies)
        thermo_group = sum(thermo_freqs_norm.get(aa, 0) for aa in group_aas)
        thermo_other = sum(thermo_freqs_norm.get(aa, 0) for aa in thermo_freqs_norm if aa not in group_aas)
        
        generic_group = sum(generic_freqs_norm.get(aa, 0) for aa in group_aas)
        generic_other = sum(generic_freqs_norm.get(aa, 0) for aa in generic_freqs_norm if aa not in group_aas)
        
        thermo_group_total += thermo_group
        thermo_other_total += thermo_other
        generic_group_total += generic_group
        generic_other_total += generic_other
    
    # Create contingency table
    contingency = np.array([
        [thermo_group_total, thermo_other_total],
        [generic_group_total, generic_other_total]
    ])
    
    # Calculate statistics
    chi2, p_value, dof, expected = chi2_contingency(contingency)
    
    # Calculate fold enrichment
    thermo_ratio = thermo_group_total / (thermo_group_total + thermo_other_total)
    generic_ratio = generic_group_total / (generic_group_total + generic_other_total)
    fold_enrichment = thermo_ratio / generic_ratio if generic_ratio > 0 else float('inf')
    
    return {
        'group': group_name,
        'thermo_freq': thermo_ratio,
        'generic_freq': generic_ratio,
        'fold_enrichment': fold_enrichment,
        'chi2': chi2,
        'p_value': p_value,
        'thermo_count': thermo_group_total,
        'generic_count': generic_group_total
    }

def main():
    if len(sys.argv) != 3:
        print("Usage: python amino_acid_enrichment.py <thermostable_file> <generic_file>")
        sys.exit(1)
    
    thermo_file = sys.argv[1]
    generic_file = sys.argv[2]
    
    # Load data
    print("Loading data...")
    thermo_data, amino_acids = load_exposed_data(thermo_file)
    generic_data, _ = load_exposed_data(generic_file)
    
    # Find common positions
    common_positions = set(thermo_data.keys()) & set(generic_data.keys())
    print(f"Analyzing {len(common_positions)} common positions")
    
    # Calculate enrichment for each group
    results = []
    for group_name, group_aas in AA_GROUPS.items():
        result = calculate_group_enrichment(thermo_data, generic_data, common_positions, group_name, group_aas)
        results.append(result)
    
    # Convert to DataFrame and sort by fold enrichment
    df_results = pd.DataFrame(results)
    df_results = df_results.sort_values('fold_enrichment', ascending=False)
    
    # Add significance indicators
    df_results['significant'] = df_results['p_value'] < 0.05
    
    # Save results
    output_file = "amino_acid_enrichment_charged_results.csv"
    df_results.to_csv(output_file, index=False)
    print(f"\nResults saved to {output_file}")
    
    # Print summary
    print(f"\nAmino acid group enrichment analysis:")
    print("=" * 80)
    print(f"{'Group':<18} {'Thermo%':<8} {'Generic%':<8} {'Fold':<8} {'P-value':<10} {'Significant'}")
    print("-" * 80)
    
    for _, row in df_results.iterrows():
        sig_mark = "***" if row['p_value'] < 0.001 else "**" if row['p_value'] < 0.01 else "*" if row['p_value'] < 0.05 else ""
        print(f"{row['group']:<18} {row['thermo_freq']:<8.3f} {row['generic_freq']:<8.3f} {row['fold_enrichment']:<8.3f} {row['p_value']:<10.6f} {sig_mark}")

if __name__ == "__main__":
    main()
