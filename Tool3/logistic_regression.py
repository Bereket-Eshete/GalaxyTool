import pandas as pd
import numpy as np
import statsmodels.api as sm
import sys

def main():
    # Parse arguments
    genotypes_path = sys.argv[1]
    phenotypes_path = sys.argv[2]
    output_tsv = sys.argv[3]
    covariates = sys.argv[4].split(',') if sys.argv[4] != 'NA' else []
    maf_thresh = float(sys.argv[5])
    missing_code = sys.argv[6]
    snp_annot_path = sys.argv[7] if len(sys.argv) > 7 else None

    # Load data
    genotypes = pd.read_csv(genotypes_path, sep='\t', index_col=0)
    phenotypes = pd.read_csv(phenotypes_path, sep='\t').set_index('sample_id')
    if snp_annot_path:
        snp_annot = pd.read_csv(snp_annot_path, sep='\t')

    # Prepare covariates matrix
    X_cov = phenotypes[covariates] if covariates else pd.DataFrame()
    if 'intercept' not in X_cov.columns:
        X_cov = sm.add_constant(X_cov)  # Add intercept

    results = []
    for snp_id in genotypes.columns:
        # Merge genotype with phenotypes/covariates
        df = phenotypes.join(genotypes[[snp_id]].rename(columns={snp_id: 'genotype'}))
        df = df[df['genotype'] != missing_code].dropna()
        
        # Check MAF
        maf = df['genotype'].mean() / 2
        if maf < maf_thresh or len(df) < 10:
            continue
        
        # Fit logistic regression
        X = pd.concat([X_cov, df['genotype']], axis=1).dropna()
        y = df.loc[X.index, 'phenotype']
        model = sm.Logit(y, X).fit(disp=0)
        
        # Store results
        results.append({
            'snp_id': snp_id,
            'beta': model.params['genotype'],
            'se': model.bse['genotype'],
            'z': model.tvalues['genotype'],
            'p_value': model.pvalues['genotype'],
            'maf': maf
        })

    # Merge with SNP annotation (if provided)
    result_df = pd.DataFrame(results)
    if snp_annot_path:
        result_df = result_df.merge(snp_annot, on='snp_id', how='left')
    
    # Save output
    result_df.to_csv(output_tsv, sep='\t', index=False)

if __name__ == "__main__":
    main()