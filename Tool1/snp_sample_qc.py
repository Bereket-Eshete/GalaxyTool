import pandas as pd
import numpy as np
import sys

def main():
    # Parse arguments from Galaxy XML
    genotypes_path = sys.argv[1]
    snp_annot_path = sys.argv[2]
    output_geno = sys.argv[3]
    output_snp = sys.argv[4]
    output_samples = sys.argv[5]
    qc_report = sys.argv[6]
    maf_thresh = float(sys.argv[7])
    snp_miss_thresh = float(sys.argv[8])
    sample_miss_thresh = float(sys.argv[9])
    missing_code = sys.argv[10]

    # Load data
    genotypes = pd.read_csv(genotypes_path, sep='\t', index_col=0)
    snp_annot = pd.read_csv(snp_annot_path, sep='\t')

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