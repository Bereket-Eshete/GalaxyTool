# Manhattan Plot Generator

This tool creates a Manhattan plot from GWAS association results and reports the top N most significant SNPs.

## Inputs

- **Association Results TSV**: A tab-separated file containing GWAS results with the following columns:

  - `snp_id`: Unique identifier for each SNP.
  - `chrom`: Chromosome (numeric or string).
  - `pos`: Position on the chromosome.
  - `p_value`: P-value for the association test.

- **P-value Significance Threshold**: Threshold for highlighting significant SNPs (default = 5e-8).

- **Number of Top Hits**: Number of top significant SNPs to report (default = 20).

## Outputs

- **Manhattan Plot PNG**: A visual representation of GWAS results, with chromosomes on the x-axis and -log10(p-value) on the y-axis.
- **Top Hits TSV**: A tab-separated file containing the top N most significant SNPs.

## Usage

Run the tool using the following command:

```bash
python manhattan_plot.py \
    gwas_data/association_results.tsv \
    Tool4/sample_outputs/ \
    5e-8 \
    20
```
