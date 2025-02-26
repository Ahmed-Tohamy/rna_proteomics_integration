# RNA-Seq and Proteomics Data Integration Script
# ------------------------------------------------------
# This script integrates RNA-Seq and proteomics data from Burkholderia thailandensis,
# computes Pearson correlation between log2 fold changes, and generates relevant visualizations.
# ------------------------------------------------------

# Required Libraries
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.preprocessing import StandardScaler

# ------------------------------------------------------
# Step 1: Load the RNA-Seq and Proteomics Data
# ------------------------------------------------------
# Expected input format:
# RNA-Seq: CSV file with columns ['GeneID', 'log2FoldChange', 'p-value']
# Proteomics: CSV file with columns ['GeneID', 'log2FoldChange', 'PValue']

rna_file = 'rna_seq_data.csv'  # Update with actual file path
proteomics_file = 'proteomics_data.csv'  # Update with actual file path

# Load data into pandas DataFrames
df_rna = pd.read_csv(rna_file)
df_pro = pd.read_csv(proteomics_file)

# Display initial data structure
print("RNA-Seq Data Preview:")
print(df_rna.head())
print("\nProteomics Data Preview:")
print(df_pro.head())

# ------------------------------------------------------
# Step 2: Preprocess and Merge Data
# ------------------------------------------------------
# Remove missing values
df_rna = df_rna.dropna()
df_pro = df_pro.dropna()

# Merge datasets on 'GeneID'
df_merged = pd.merge(df_rna, df_pro, on='GeneID', suffixes=('_RNA', '_Pro'))

# Display merged data
print("\nMerged Data Preview:")
print(df_merged.head())

# ------------------------------------------------------
# Step 3: Normalize Data Using StandardScaler
# ------------------------------------------------------
scaler = StandardScaler()
df_merged['log2FoldChange_RNA_norm'] = scaler.fit_transform(df_merged[['log2FoldChange_RNA']])
df_merged['log2FoldChange_Pro_norm'] = scaler.fit_transform(df_merged[['log2FoldChange_Pro']])

# ------------------------------------------------------
# Step 4: Compute Pearson Correlation
# ------------------------------------------------------
correlation, p_value = stats.pearsonr(df_merged['log2FoldChange_RNA_norm'], df_merged['log2FoldChange_Pro_norm'])
print(f"\nPearson Correlation Coefficient: {correlation:.2f}")
print(f"P-Value: {p_value:.2e}")

# ------------------------------------------------------
# Step 5: Generate Visualizations
# ------------------------------------------------------

# Scatter Plot with Regression Line
plt.figure(figsize=(10, 6))
sns.regplot(x='log2FoldChange_RNA_norm', y='log2FoldChange_Pro_norm', data=df_merged, scatter_kws={'s': 30}, line_kws={'color': 'red'})
plt.xlabel('Normalized Log2 Fold Change (RNA-Seq)')
plt.ylabel('Normalized Log2 Fold Change (Proteomics)')
plt.title('Correlation Between RNA-Seq and Proteomics Data')
plt.grid(True)
plt.savefig('scatter_plot.png', dpi=300)
plt.show()

# Density Plot
plt.figure(figsize=(10, 6))
sns.kdeplot(df_merged['log2FoldChange_RNA_norm'], label='RNA-Seq', fill=True)
sns.kdeplot(df_merged['log2FoldChange_Pro_norm'], label='Proteomics', fill=True)
plt.xlabel('Normalized Log2 Fold Change')
plt.ylabel('Density')
plt.title('Distribution of Normalized Log2 Fold Changes')
plt.legend()
plt.savefig('density_plot.png', dpi=300)
plt.show()

# Violin Plot
plt.figure(figsize=(10, 6))
sns.violinplot(data=[df_merged['log2FoldChange_RNA_norm'], df_merged['log2FoldChange_Pro_norm']], inner='quartile')
plt.xticks([0, 1], ['RNA-Seq', 'Proteomics'])
plt.ylabel('Normalized Log2 Fold Change')
plt.title('Violin Plot of Normalized Log2 Fold Changes')
plt.savefig('violin_plot.png', dpi=300)
plt.show()

# ------------------------------------------------------
# Step 6: Save Processed Data
# ------------------------------------------------------
df_merged.to_csv('integrated_data.csv', index=False)
print("\nProcessed data saved as 'integrated_data.csv'")

# ------------------------------------------------------
# Notes:
# - This script requires Pandas, NumPy, SciPy, scikit-learn, Seaborn, and Matplotlib.
# - To install dependencies, run: pip install pandas numpy scipy scikit-learn seaborn matplotlib
# - Ensure that the RNA-Seq and proteomics input files are formatted correctly before execution.
# ------------------------------------------------------
