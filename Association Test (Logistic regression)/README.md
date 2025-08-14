# Logistic Regression Tool

This tool performs logistic regression for each SNP in a genotype dataset, adjusting for optional covariates.

## Inputs
- **Genotypes TSV**: A tab-separated file with samples as rows and SNPs as columns (values: 0/1/2).
- **Phenotypes TSV**: A tab-separated file with `sample_id` and `phenotype` columns.
- **Covariates**: Optional comma-separated list of covariate column names.
- **SNP Annotation TSV**: A tab-separated file with SNP metadata (`snp_id`, `chrom`, `pos`, `maf`).

## Parameters
- **MAF Threshold**: Minimum minor allele frequency (default = 0.01).
- **Missing Value Code**: Code used to represent missing values (default = "NA").

## Outputs
- **Logistic Regression Results TSV**: A tab-separated file containing:
  - `snp_id`: Unique identifier for each SNP.
  - `chrom`: Chromosome.
  - `pos`: Position.
  - `beta`: Effect size (log odds ratio).
  - `se`: Standard error of the effect size.
  - `z`: Z-score.
  - `p_value`: P-value.
  - `maf`: Minor allele frequency.

## Usage
Run the tool using the following command:
```bash
python logistic_regression.py \
    gwas_data/genotypes.tsv \
    gwas_data/phenotypes.tsv \
    Tool3/sample_outputs/ \
    age,sex \
    0.01 \
    NA \
    gwas_data/snp_annotation.tsv