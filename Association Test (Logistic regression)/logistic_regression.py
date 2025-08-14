#!/usr/bin/env python3
import pandas as pd
import numpy as np
import statsmodels.api as sm
import sys
import os

def main():
    # Parse arguments
    genotypes_path = sys.argv[1]
    phenotypes_path = sys.argv[2]
    output_dir = sys.argv[3]  # Directory for outputs
    covariates = sys.argv[4].split(',') if sys.argv[4] != 'NA' else []
    maf_thresh = float(sys.argv[5])
    missing_code = sys.argv[6]
    snp_annot_path = sys.argv[7]  # Now REQUIRED

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Define output path
    output_tsv = os.path.join(output_dir, "logistic_regression_results.tsv")

    # Load data
    genotypes = pd.read_csv(genotypes_path, sep='\t', index_col=0)
    phenotypes = pd.read_csv(phenotypes_path, sep='\t').set_index('sample_id')
    snp_annot = pd.read_csv(snp_annot_path, sep='\t')

    # Verify SNP annotation has required columns
    required_columns = ['snp_id', 'chrom', 'pos', 'maf']
    if not all(col in snp_annot.columns for col in required_columns):
        raise ValueError(f"SNP annotation must contain: {required_columns}")

    # Prepare covariates matrix
    X_cov = phenotypes[covariates] if covariates else pd.DataFrame()
    if 'intercept' not in X_cov.columns:
        X_cov = sm.add_constant(X_cov)  # Add intercept

    results = []
    for snp_id in genotypes.columns:
        # Skip if SNP not in annotation
        if snp_id not in snp_annot['snp_id'].values:
            print(f"Warning: SNP {snp_id} not in annotation. Skipping.", file=sys.stderr)
            continue

        # Get chrom/pos/maf from annotation
        snp_data = snp_annot.loc[snp_annot['snp_id'] == snp_id, ['chrom', 'pos', 'maf']].iloc[0]
        chrom, pos, annot_maf = snp_data['chrom'], snp_data['pos'], snp_data['maf']

        # Skip if MAF below threshold
        if annot_maf < maf_thresh:
            continue

        # Merge genotype with phenotypes
        df = phenotypes.join(genotypes[[snp_id]].rename(columns={snp_id: 'genotype'}))
        df = df[df['genotype'] != missing_code].dropna()

        # Skip if <10 samples
        if len(df) < 10:
            continue

        # Fit logistic regression
        X = pd.concat([X_cov, df['genotype']], axis=1).dropna()
        y = df.loc[X.index, 'phenotype']
        model = sm.Logit(y, X).fit(disp=0)

        results.append({
            'snp_id': snp_id,
            'chrom': chrom,
            'pos': pos,
            'beta': model.params['genotype'],
            'se': model.bse['genotype'],
            'z': model.tvalues['genotype'],
            'p_value': model.pvalues['genotype'],
            'maf': annot_maf
        })

    # Save output
    result_df = pd.DataFrame(results)
    result_df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Output saved to {output_tsv}")

if __name__ == "__main__":
    main()