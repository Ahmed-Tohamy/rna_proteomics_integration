# **RNA-Seq and Proteomics Data Integration**

This repository contains a Python script for integrating **RNA-Seq** and **proteomics** data from **Burkholderia thailandensis**, computing **Pearson correlation** between log2 fold changes, and generating various **visualizations**.

---

##  Features
- **Data Integration**: Merges RNA-Seq and proteomics datasets based on `GeneID`.
- **Data Preprocessing**: Removes missing values and normalizes data using **StandardScaler**.
- **Correlation Analysis**: Computes **Pearson correlation** between RNA-Seq and proteomics log2 fold changes.
- **Visualization**: Generates:
  - **Scatter plot with regression line** (Correlation analysis)
  - **Density plot** (Distribution of normalized log2 fold changes)
  - **Violin plot** (Comparison of RNA-Seq and Proteomics fold changes)
- **Processed Data Storage**: Saves the integrated and normalized dataset as `integrated_data.csv`.

---

## Input Data Format
This script expects **two CSV files** with the following columns:

1. **RNA-Seq Data (`rna_seq_data.csv`)**
GeneID, log2FoldChange, p-value

2. **Proteomics Data (`proteomics_data.csv`)**
GeneID, log2FoldChange, PValue


---

## Installation
To use this script, install the required dependencies:

```sh
pip install pandas numpy scipy scikit-learn seaborn matplotlib


To execute the script, use:
python RNA-Seq_Proteomics_Integration.py
