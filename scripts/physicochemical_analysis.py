#!/usr/bin/env python3
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, mannwhitneyu
import sys

### Calculate the values of different phyiscochemical properties in every common position of two sets of PSSMs using the PDB position identifiers. Run MW-U and T-test.

# Physicochemical property scales
PROPERTIES = {
    'Hydrophobicity_KD': {  # Kyte-Doolittle scale
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    },
    'Volume': {  # Amino acid volume (Å³)
        'A': 88.6, 'R': 173.4, 'N': 114.1, 'D': 111.1, 'C': 108.5,
        'Q': 143.8, 'E': 138.4, 'G': 60.1, 'H': 153.2, 'I': 166.7,
        'L': 166.7, 'K': 168.6, 'M': 162.9, 'F': 189.9, 'P': 112.7,
        'S': 89.0, 'T': 116.1, 'W': 227.8, 'Y': 193.6, 'V': 140.0
    },
    'Charge': {  # Formal charge at pH 7 (IMGT-consistent)
        'A': 0, 'R': 1, 'N': 0, 'D': -1, 'C': 0,
        'Q': 0, 'E': -1, 'G': 0, 'H': 1, 'I': 0,
        'L': 0, 'K': 1, 'M': 0, 'F': 0, 'P': 0,
        'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0
    },
    'Polarity': {  # Zimmerman polarity scale
        'A': 0.0, 'R': 52.0, 'N': 3.38, 'D': 49.7, 'C': 1.48,
        'Q': 3.53, 'E': 49.9, 'G': 0.0, 'H': 51.6, 'I': 0.13,
        'L': 0.13, 'K': 49.5, 'M': 1.43, 'F': 0.35, 'P': 1.58,
        'S': 1.67, 'T': 1.66, 'W': 2.1, 'Y': 1.61, 'V': 0.13
    },
    'Flexibility': {  # Bhaskaran & Ponnuswamy (ProtScale)
        'A': 0.360, 'R': 0.530, 'N': 0.460, 'D': 0.510, 'C': 0.350,
        'Q': 0.490, 'E': 0.500, 'G': 0.540, 'H': 0.320, 'I': 0.460,
        'L': 0.370, 'K': 0.470, 'M': 0.300, 'F': 0.310, 'P': 0.510,
        'S': 0.510, 'T': 0.440, 'W': 0.310, 'Y': 0.420, 'V': 0.390
    }
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

def calculate_weighted_property(freq_dict, property_scale):
    """Calculate weighted average property value"""
    total = 0
    weight_sum = 0
    for aa, freq in freq_dict.items():
        total += freq * property_scale[aa]
        weight_sum += freq
    return total / weight_sum if weight_sum > 0 else 0

def main():
    if len(sys.argv) != 3:
        print("Usage: python physicochemical_analysis.py <thermostable_file> <generic_file>")
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
    
    # Calculate properties for each position
    results = []
    for pos in sorted(common_positions):
        result = {'PDB_resnum': pos}
        
        for prop_name, prop_scale in PROPERTIES.items():
            thermo_value = calculate_weighted_property(thermo_data[pos], prop_scale)
            generic_value = calculate_weighted_property(generic_data[pos], prop_scale)
            
            result[f'thermo_{prop_name}'] = thermo_value
            result[f'generic_{prop_name}'] = generic_value
            result[f'diff_{prop_name}'] = thermo_value - generic_value
        
        results.append(result)
    
    df_results = pd.DataFrame(results)
    
    # Calculate overall statistics
    summary = []
    for prop_name in PROPERTIES.keys():
        thermo_col = f'thermo_{prop_name}'
        generic_col = f'generic_{prop_name}'
        diff_col = f'diff_{prop_name}'
        
        thermo_vals = df_results[thermo_col]
        generic_vals = df_results[generic_col]
        diff_vals = df_results[diff_col]
        
        # Statistical tests
        t_stat, t_pvalue = ttest_ind(thermo_vals, generic_vals)
        u_stat, u_pvalue = mannwhitneyu(thermo_vals, generic_vals, alternative='two-sided')
        
        summary.append({
            'Property': prop_name,
            'Thermo_mean': thermo_vals.mean(),
            'Generic_mean': generic_vals.mean(),
            'Difference_mean': diff_vals.mean(),
            'Difference_std': diff_vals.std(),
            't_test_pvalue': t_pvalue,
            'mannwhitney_pvalue': u_pvalue
        })
    
    df_summary = pd.DataFrame(summary)
    
    # Save results
    df_results.to_csv("physicochemical_detailed.csv", index=False)
    df_summary.to_csv("physicochemical_summary.csv", index=False)
    
    print(f"\nPhysicochemical property comparison:")
    print("=" * 80)
    print(f"{'Property':<15} {'Thermo':<8} {'Generic':<8} {'Diff':<8} {'t-test p':<10} {'MW p':<10}")
    print("-" * 80)
    
    for _, row in df_summary.iterrows():
        t_sig = "***" if row['t_test_pvalue'] < 0.001 else "**" if row['t_test_pvalue'] < 0.01 else "*" if row['t_test_pvalue'] < 0.05 else ""
        mw_sig = "***" if row['mannwhitney_pvalue'] < 0.001 else "**" if row['mannwhitney_pvalue'] < 0.01 else "*" if row['mannwhitney_pvalue'] < 0.05 else ""
        
        print(f"{row['Property']:<15} {row['Thermo_mean']:<8.3f} {row['Generic_mean']:<8.3f} {row['Difference_mean']:<8.3f} {row['t_test_pvalue']:<10.6f}{t_sig} {row['mannwhitney_pvalue']:<10.6f}{mw_sig}")

if __name__ == "__main__":
    main()
