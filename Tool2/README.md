# Allele Frequency Calculator

This tool calculates allele counts and minor allele frequencies (MAF) for each SNP in a genotype dataset.

## Inputs

- **Genotypes TSV**: A tab-separated file with samples as rows and SNPs as columns (values: 0/1/2).
- **Optional SNP Annotation TSV**: A tab-separated file with SNP metadata (`snp_id`, `chrom`, `pos`).
- **Optional Sample Subset TSV**: A tab-separated file containing a list of sample IDs to include in the analysis.

## Parameters

- **Missing Value Code**: Code used to represent missing values (default: "NA").

## Outputs

- **Allele Frequencies TSV**: A tab-separated file containing:
  - `snp_id`: Unique identifier for each SNP.
  - `count_0`, `count_1`, `count_2`: Number of samples with each allele count.
  - `n_non_missing`: Number of non-missing samples for the SNP.
  - `MAF`: Minor allele frequency.

## Usage

Run the tool using the following command:

```bash
python allele_frequency_calculator.py \
    gwas_data/genotypes.tsv \
    Tool2/sample_outputs/ \
    NA \
    gwas_data/sample_subset.tsv \
    gwas_data/snp_annotation.tsv
```
