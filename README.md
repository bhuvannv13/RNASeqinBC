# RNASeqinBC
analyse a publicly available breast cancer dataset. Breast cancer is the most diagnosed cancer in women. There are several subtypes of diseases characterised by different genetic aberrations drivers. The Human Epidermal growth factor Receptor 2 (HER2) amplified, also termed as ERBB2+, breast cancer is one of the most aggressive subtypes. 


---

# RNASeq Data Analysis: HER2 Amplification

This repository contains the analysis of RNA sequencing data to study HER2 amplification's impact on gene expression and survival in cancer patients.

## Table of Contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Key Steps](#key-steps)
5. [Results](#results)
6. [Dependencies](#dependencies)
7. [License](#license)

## Overview

This project identifies differentially expressed genes, performs pathway enrichment, and examines survival using RNA-Seq data. The analysis focuses on HER2 amplification in cancer samples.

## Installation

1. Dataset:
  https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018

3. Install R and RStudio if not already installed.

4. Ensure the required R packages are available (see [Dependencies](#dependencies)).

## Usage

1. Set the paths to the input data files:
   - `data_clinical_patient.txt`: Clinical data.
   - `data_mrna_seq_v2_rsem.txt`: RNA-Seq data.
   - `data_cna.txt`: Copy number alteration data.

2. Modify the `folder_path` variable in the notebook to point to the directory containing the above files.

3. Open the `.Rmd` file in RStudio and run all cells to execute the analysis.

## Key Steps

1. **Preprocessing**:
   - Load clinical, RNA-Seq, and CNA datasets.
   - Map patient IDs and generate metadata for HER2 amplification status and cancer stage.

2. **Normalization and Differential Expression**:
   - Use `DESeq2` to normalize RNA-Seq data.
   - Identify differentially expressed genes based on HER2 amplification.

3. **Pathway Enrichment**:
   - Perform Gene Ontology enrichment analysis on over-expressed genes.

4. **Visualization**:
   - PCA for global expression trends.
   - Heatmap for differentially expressed genes.

5. **Survival Analysis**:
   - Apply Lasso regularization with Cox regression to identify survival predictors.

## Results

- **Differential Expression**:
  Top genes associated with HER2 amplification.
- **Pathway Enrichment**:
  Biological processes enriched among over-expressed genes.
- **Survival Predictors**:
  Key genes impacting overall survival.

## Dependencies

The following R packages are required:

- `DESeq2`
- `clusterProfiler`
- `org.Hs.eg.db`
- `pheatmap`
- `glmnet`
- `survival` (for Cox regression)

To install them:
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db"))
install.packages(c("pheatmap", "glmnet", "survival"))
```



Feel free to update the repository link and add more details as needed. Let me know if you'd like this drafted into a file!
