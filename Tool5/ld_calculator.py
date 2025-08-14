#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from scipy.stats import pearsonr

def calculate_ld(genotypes):
    """Calculate pairwise r² for all SNPs"""
    n_snps = genotypes.shape[1]
    ld_matrix = np.full((n_snps, n_snps), np.nan)  # Initialize with NaN
    
    for i in range(n_snps):
        for j in range(i, n_snps):
            # Skip constant values (monomorphic SNPs)
            if genotypes.iloc[:, i].nunique() == 1 or genotypes.iloc[:, j].nunique() == 1:
                r2 = np.nan
            else:
                try:
                    r, _ = pearsonr(genotypes.iloc[:, i], genotypes.iloc[:, j])
                    r2 = r ** 2
                except ValueError:  # Handle cases where Pearson correlation fails
                    r2 = np.nan
            ld_matrix[i, j] = ld_matrix[j, i] = r2
            
    return pd.DataFrame(ld_matrix, 
                        index=genotypes.columns, 
                        columns=genotypes.columns)

def main():
    # Verify arguments
    if len(sys.argv) != 9:
        print("Usage: python ld_calculator.py <genotypes> <snp_annot> <focal_snp> <window_kb> <min_maf> <missing_code> <output_matrix> <output_heatmap>")
        sys.exit(1)

    try:
        # Load data
        genotypes = pd.read_csv(sys.argv[1], sep='\t', index_col=0, dtype=str)
        
        # Load SNP annotation file
        snp_annot = pd.read_csv(sys.argv[2], sep='\t')
        snp_annot.columns = snp_annot.columns.str.strip()  # Clean column names
        snp_annot = snp_annot.astype({'snp_id': str, 'chrom': str, 'pos': int, 'maf': float})  # Set dtypes
        
        # Debugging output
        print("Column names in snp_annotation.tsv:", snp_annot.columns.tolist())
        
        # Parse arguments
        focal_snp = sys.argv[3]
        window_kb = float(sys.argv[4])
        min_maf = float(sys.argv[5])
        missing_code = sys.argv[6]
        output_matrix = sys.argv[7]
        output_heatmap = sys.argv[8]

        # Get focal SNP position
        focal_match = snp_annot[snp_annot['snp_id'] == focal_snp]
        if focal_match.empty:
            raise ValueError(f"Focal SNP {focal_snp} not found in annotation")
        print("Focal SNP match:", focal_match)
        chrom, center_pos = focal_match[['chrom', 'pos']].values[0]
        print(f"Processing SNP: {focal_snp} at {chrom}:{center_pos}")

        # Select SNPs in window
        window_bp = int(window_kb * 1000)
        window_snps = snp_annot[
            (snp_annot['chrom'] == chrom) & 
            (snp_annot['pos'].between(center_pos - window_bp // 2, center_pos + window_bp // 2))
        ]['snp_id']
        
        if window_snps.empty:
            raise ValueError(f"No SNPs in {window_kb}kb window around {focal_snp}")

        # Filter genotypes
        available_snps = list(set(window_snps) & set(genotypes.columns))
        if len(available_snps) < 2:
            raise ValueError(f"Only {len(available_snps)} SNP(s) available after merging annotation and genotypes")
        
        geno_window = genotypes[available_snps].replace(missing_code, np.nan)
        geno_window = geno_window.apply(pd.to_numeric, errors='coerce')  # Convert to numeric, coercing errors
        geno_window = geno_window.dropna(how='all').dropna(axis=1, how='all')
        
        # MAF filtering
        mafs = geno_window.mean(axis=0) / 2
        geno_window = geno_window.loc[:, mafs >= min_maf]
        
        if geno_window.shape[1] < 2:
            raise ValueError(
                f"Only {geno_window.shape[1]} SNP(s) passed MAF filter (≥{min_maf})\n"
                f"Try: (1) Lower MAF threshold (2) Larger window size"
            )

        # Calculate LD matrix
        print(f"Calculating LD for {geno_window.shape[1]} SNPs...")
        ld_matrix = calculate_ld(geno_window)
        
        # Save outputs
        ld_matrix.to_csv(output_matrix, sep='\t', float_format='%.4f')
        print(f"Saved LD matrix to {output_matrix}")

        # Create heatmap
        plt.figure(figsize=(10, 8))
        sns.heatmap(
            ld_matrix,
            cmap='coolwarm',
            vmin=0,
            vmax=1,
            square=True,
            xticklabels=ld_matrix.columns,
            yticklabels=ld_matrix.index
        )
        plt.title(f"LD around {focal_snp} ({window_kb}kb window)\n{geno_window.shape[1]} SNPs, MAF ≥ {min_maf}")
        plt.tight_layout()
        plt.savefig(output_heatmap, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved heatmap to {output_heatmap}")

    except Exception as e:
        print(f"\nERROR: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()