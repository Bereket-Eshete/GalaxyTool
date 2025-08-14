GWAS Tools Suite
This repository contains a collection of tools designed for Genome-Wide Association Studies (GWAS). Each tool performs a specific task in the GWAS pipeline, such as quality control, association analysis, visualization, and more. These tools are modular, well-documented, and can be used independently or as part of a larger workflow.

Overview of Tools
The suite includes the following tools:

Tool
1: SNP and Sample QC
Performs quality control on genotype data by filtering SNPs and samples based on missingness and minor allele frequency (MAF) thresholds.
Folder: Tool1/
Tool 2: Allele Frequency Calculator
Computes per-SNP allele counts and minor allele frequencies (MAF) from genotype data.
Folder: Tool2/
Tool 3: Logistic Regression
Performs logistic regression for each SNP, adjusting for optional covariates, to identify associations with phenotypes.
Folder: Tool3/
Tool 4: Manhattan Plot Generator
Creates a Manhattan plot from GWAS results and reports the top N most significant SNPs.
Folder: Tool4/
Tool 5: Windowed LD Calculator
Calculates pairwise linkage disequilibrium (LD) for SNPs within a specified genomic window around a focal SNP.
Folder: Tool5/
Each tool has its own folder containing:

Python script (*.py)
XML wrapper for Galaxy integration (*.xml)
Documentation (README.md)
Example outputs stored in the sample_outputs/ folder
Getting Started
Prerequisites
Python 3.x
Required libraries: pandas, numpy, matplotlib, seaborn, statsmodels, scipy
Input files in tab-separated format (TSV)
Installation
1.Clone this repository:
bash
git clone https://github.com/your-repo/gwas-tools.git
cd gwas-tools
2.Install dependencies:
bash
pip install pandas numpy matplotlib seaborn statsmodels scipy
Explore each tool's folder for detailed instructions and usage examples.
Usage
To use a specific tool, navigate to its folder and follow the instructions in the README.md file. For example:
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
Folder Structure 
gwas-tools/
├── SNP and Sample QC/              
├── Allele Frequency Calculator/           
├── Logistic Regression/             
├── Manhattan Plot Generator/             
├── Windowed LD Calculator          
├── README.md          
└── gwas_data/      
