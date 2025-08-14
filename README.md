# GWAS Tools Suite

This repository contains a collection of tools designed for **Genome-Wide Association Studies (GWAS)**.  
Each tool performs a specific task in the GWAS pipeline, such as **quality control**, **association analysis**, **visualization**, and more.  
These tools are **modular**, **well-documented**, and can be used **independently** or as part of a larger workflow.

---

## 📌 Overview of Tools

| Tool | Description | Folder |
|------|-------------|--------|
| **Tool 1: SNP and Sample QC** | Performs quality control on genotype data by filtering SNPs and samples based on missingness and minor allele frequency (MAF) thresholds. | `Tool1/` |
| **Tool 2: Allele Frequency Calculator** | Computes per-SNP allele counts and minor allele frequencies (MAF) from genotype data. | `Tool2/` |
| **Tool 3: Logistic Regression** | Performs logistic regression for each SNP, adjusting for optional covariates, to identify associations with phenotypes. | `Tool3/` |
| **Tool 4: Manhattan Plot Generator** | Creates a Manhattan plot from GWAS results and reports the top N most significant SNPs. | `Tool4/` |
| **Tool 5: Windowed LD Calculator** | Calculates pairwise linkage disequilibrium (LD) for SNPs within a specified genomic window around a focal SNP. | `Tool5/` |

---

## 📂 Folder Structure

gwas-tools/
├── SNP and Sample QC/
├── Allele Frequency Calculator/
├── Logistic Regression/
├── Manhattan Plot Generator/
├── Windowed LD Calculator/
├── README.md
└── gwas_data/

---

## 🚀 Getting Started

### **Prerequisites**
- **Python 3.x**
- Required Python libraries:
  - `pandas`
  - `numpy`
  - `matplotlib`
  - `seaborn`
  - `statsmodels`
  - `scipy`
- Input files in **tab-separated format (TSV)**

---

### **Installation**

```bash
# 1. Clone this repository
git clone https://github.com/your-repo/gwas-tools.git
cd gwas-tools

# 2. Install dependencies
pip install pandas numpy matplotlib seaborn statsmodels scipy
Usage

To use a specific tool, navigate to its folder and follow the instructions in the README.md file for that tool.

Example: Running Tool 1 (SNP and Sample QC)
cd Tool1/
python snp_sample_qc.py \
    gwas_data/genotypes.tsv \
    gwas_data/snp_annotation.tsv \
    Tool1/sample_outputs/ \
    0.01 \
    0.05 \
    0.1 \
    NA
All outputs will be saved in the respective sample_outputs/ folder.
