# SNP and Sample QC Tool

This tool performs quality control (QC) on genotype data by filtering SNPs and samples based on missingness and minor allele frequency (MAF) thresholds.

## Inputs

- **Genotypes TSV**: A tab-separated file with samples as rows and SNPs as columns (values: 0/1/2).
- **SNP Annotation TSV**: A tab-separated file with SNP metadata (`snp_id`, `chrom`, `pos`, `maf`).

## Parameters

- **MAF Threshold**: Minimum minor allele frequency (default: 0.01).
- **SNP Missingness Threshold**: Maximum proportion of missing values allowed for SNPs (default: 0.05).
- **Sample Missingness Threshold**: Maximum proportion of missing values allowed for samples (default: 0.1).
- **Missing Value Code**: Code used to represent missing values (default: "NA").

## Outputs

- **Filtered Genotypes**: Genotype data after QC.
- **Filtered SNP Annotation**: SNP metadata after QC.
- **Filtered Samples**: List of retained samples.
- **QC Report**: Summary of QC metrics.

## Usage

Run the tool using the following command:

```bash
python snp_sample_qc.py \
    genotypes.tsv \
    snp_annotation.tsv \
    sample_outputs/ \
    0.01 \
    0.05 \
    0.1 \
    NA
```
