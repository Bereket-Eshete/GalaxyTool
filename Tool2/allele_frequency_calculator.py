import pandas as pd
import numpy as np
import sys

def main():
    # Parse arguments
    genotypes_path = sys.argv[1]
    output_tsv = sys.argv[2]
    missing_code = sys.argv[3]  # e.g., "NA" or "-9"
    sample_subset_path = sys.argv[4] if len(sys.argv) > 4 else None
    snp_annot_path = sys.argv[5] if len(sys.argv) > 5 else None

    # Load data
    genotypes = pd.read_csv(genotypes_path, sep='\t', index_col=0)
    if sample_subset_path:
        sample_subset = pd.read_csv(sample_subset_path, sep='\t')['sample_id'].tolist()
        genotypes = genotypes.loc[sample_subset]

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
    if snp_annot_path:
        snp_annot = pd.read_csv(snp_annot_path, sep='\t')
        result_df = result_df.merge(snp_annot, on='snp_id', how='left')
    
    # Save output
    result_df.to_csv(output_tsv, sep='\t', index=False)

if __name__ == "__main__":
    main()