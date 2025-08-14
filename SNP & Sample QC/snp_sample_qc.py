#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys
import os

def main():
    # Parse arguments from Galaxy XML
    genotypes_path = sys.argv[1]
    snp_annot_path = sys.argv[2]
    output_dir = sys.argv[3]  # Directory for outputs
    maf_thresh = float(sys.argv[4])
    snp_miss_thresh = float(sys.argv[5])
    sample_miss_thresh = float(sys.argv[6])
    missing_code = sys.argv[7]

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Define output paths
    output_geno = os.path.join(output_dir, "filtered_genotypes.tsv")
    output_snp = os.path.join(output_dir, "filtered_snp_annotation.tsv")
    output_samples = os.path.join(output_dir, "filtered_samples.tsv")
    qc_report = os.path.join(output_dir, "qc_report.txt")

    # Load genotype data
    genotypes = pd.read_csv(genotypes_path, sep='\t', index_col=0)

    # Load SNP annotation file
    snp_annot = pd.read_csv(snp_annot_path, sep='\t')
    snp_annot.columns = snp_annot.columns.str.strip()  # Remove leading/trailing whitespace
    
    # Debugging output
    print("Cleaned column names:", snp_annot.columns.tolist())
    
    # Validate required columns
    required_columns = ['snp_id', 'chrom', 'pos', 'maf']
    if not all(col in snp_annot.columns for col in required_columns):
        raise ValueError(f"snp_annotation.tsv is missing required columns: {required_columns}")
    
    # Ensure proper data types
    try:
        snp_annot['maf'] = snp_annot['maf'].astype(float)
    except ValueError:
        raise ValueError("The 'maf' column contains non-numeric values. Please ensure all values are numeric.")

    # QC steps
    # 1. Filter SNPs by missingness
    snp_missingness = genotypes.apply(lambda x: (x == missing_code).mean(), axis=0)
    snp_pass_miss = snp_missingness[snp_missingness <= snp_miss_thresh].index
    genotypes = genotypes[snp_pass_miss]

    # 2. Filter SNPs by MAF (from annotation or recalculate)
    snp_annot_filtered = snp_annot[snp_annot['snp_id'].isin(snp_pass_miss)]
    snp_annot_filtered = snp_annot_filtered[snp_annot_filtered['maf'] >= maf_thresh]
    genotypes = genotypes[snp_annot_filtered['snp_id']]

    # 3. Filter samples by missingness
    sample_missingness = genotypes.apply(lambda x: (x == missing_code).mean(), axis=1)
    sample_pass_miss = sample_missingness[sample_missingness <= sample_miss_thresh].index
    genotypes = genotypes.loc[sample_pass_miss]

    # Save outputs
    genotypes.to_csv(output_geno, sep='\t')
    snp_annot_filtered.to_csv(output_snp, sep='\t', index=False)
    pd.DataFrame({'sample_id': sample_pass_miss}).to_csv(output_samples, sep='\t', index=False)

    # Generate QC report
    with open(qc_report, 'w') as f:
        f.write(f"Initial SNPs: {len(snp_missingness)}\n")
        f.write(f"SNPs after missingness filter: {len(snp_pass_miss)}\n")
        f.write(f"SNPs after MAF filter: {len(snp_annot_filtered)}\n")
        f.write(f"Samples retained: {len(sample_pass_miss)}\n")

if __name__ == "__main__":
    main()