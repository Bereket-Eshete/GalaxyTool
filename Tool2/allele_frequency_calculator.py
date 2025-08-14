#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
import os

def main():
    # Parse arguments
    genotypes_path = sys.argv[1]
    output_dir = sys.argv[2]  # Directory for outputs
    missing_code = sys.argv[3]  # e.g., "NA" or "-9"
    sample_subset_path = sys.argv[4] if len(sys.argv) > 4 else None
    snp_annot_path = sys.argv[5] if len(sys.argv) > 5 else None

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Define output path
    output_tsv = os.path.join(output_dir, "allele_frequencies.tsv")

    # Load genotype data
    genotypes = pd.read_csv(genotypes_path, sep='\t', index_col=0)

    # Load sample subset (if provided)
    if sample_subset_path and os.path.exists(sample_subset_path):
        try:
            sample_subset = pd.read_csv(sample_subset_path, sep='\t')['sample_id'].tolist()
            genotypes = genotypes.loc[sample_subset]
            print(f"Filtered genotypes using sample subset from {sample_subset_path}")
        except Exception as e:
            print(f"Warning: Could not load sample subset file '{sample_subset_path}'. Error: {e}")
    else:
        print("No sample subset file provided or file does not exist. Using all samples.")

    # Calculate allele counts
    results = []
    for snp_id in genotypes.columns:
        snp_data = genotypes[snp_id]
        snp_data = snp_data[snp_data != missing_code]  # Drop missing
        n_non_missing = len(snp_data)
        
        if n_non_missing < 10:
            counts = {'count_0': np.nan, 'count_1': np.nan, 'count_2': np.nan}
            maf = np.nan
        else:
            counts = {
                'count_0': (snp_data == 0).sum(),
                'count_1': (snp_data == 1).sum(),
                'count_2': (snp_data == 2).sum()
            }
            maf = (counts['count_1'] + 2 * counts['count_2']) / (2 * n_non_missing)
        
        results.append({
            'snp_id': snp_id,
            'count_0': counts['count_0'],
            'count_1': counts['count_1'],
            'count_2': counts['count_2'],
            'n_non_missing': n_non_missing,
            'MAF': maf
        })

    # Merge with SNP annotation (if provided)
    result_df = pd.DataFrame(results)
    if snp_annot_path and os.path.exists(snp_annot_path):
        try:
            snp_annot = pd.read_csv(snp_annot_path, sep='\t')
            result_df = result_df.merge(snp_annot, on='snp_id', how='left')
            print(f"Merged SNP annotation from {snp_annot_path}")
        except Exception as e:
            print(f"Warning: Could not load SNP annotation file '{snp_annot_path}'. Error: {e}")
    else:
        print("No SNP annotation file provided or file does not exist. Skipping merge.")

    # Save output
    result_df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Output saved to {output_tsv}")

if __name__ == "__main__":
    main()