"""
Generate a boxplot comparing Karlin δ* values across biogenic amine (BA) gene types,
overlay individual data points, and annotate statistically significant pairwise
differences based on precomputed Wilcoxon rank-sum tests with FDR correction.

Input:
- 'karlin_delta_star_results.csv' containing 'Gene_Name' and 'Delta_Star'
- 'ba_delta_star_pairwise_stats.csv' containing pairwise test results with adjusted p-values

Output:
- High-resolution PNG plot saved as 'ba_gene_delta_star_boxplot.png'

Usage:
- Ensure the input CSV files are present in the working directory.
- Run this script to produce the plot visualizing compositional bias differences.

Requirements:
- pandas
- seaborn
- matplotlib
- numpy
Author:
    Aqib Javaid
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Load data
df = pd.read_csv("karlin_delta_star_results.csv")
comparison_df = pd.read_csv("ba_delta_star_pairwise_stats.csv")

plt.figure(figsize=(10, 6))
sns.set(style="whitegrid")

# Define color palette
palette = sns.color_palette("deep")

# Boxplot
sns.boxplot(
    x='Gene_Name',
    y='Delta_Star',
    data=df,
    palette=palette,
    showfliers=False
)

# Overlay jittered points
sns.stripplot(
    x='Gene_Name',
    y='Delta_Star',
    data=df,
    color='#444444',  # Dark gray
    alpha=0.6,
    jitter=True,
    size=5
)

# Significance bars
gene_order = df['Gene_Name'].unique()
y_max = df['Delta_Star'].max()
height = y_max * 0.05
offset = y_max * 0.02

for i, row in comparison_df.iterrows():
    if not row['Significant']:
        continue
    g1, g2 = row['Gene1'], row['Gene2']
    x1 = np.where(gene_order == g1)[0][0]
    x2 = np.where(gene_order == g2)[0][0]
    x1, x2 = min(x1, x2), max(x1, x2)
    y = y_max + offset + i * height

    plt.plot([x1, x1, x2, x2], [y, y + height/2, y + height/2, y], lw=1.5, c='black')

    p = row['p_adj']
    if p < 0.001:
        stars = '***'
    elif p < 0.01:
        stars = '**'
    elif p < 0.05:
        stars = '*'
    else:
        stars = 'ns'

    # Custom star text
    plt.text(
        (x1 + x2) / 2,
        y + height/2 + 0.001,  # offset above bar
        stars,
        ha='center',
        va='bottom',
        fontsize=16,
        color='black',
        fontweight='bold'
    )

# Finalize and save
plt.xlabel("BA Gene", fontsize=13, weight='bold')
plt.ylabel(r"δ* (Karlin Dinucleotide Signature Difference)", fontsize=13, weight='bold')
plt.title("Compositional Bias of BA Genes Across LAB Genomes", fontsize=14, weight='bold')
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
plt.tight_layout()
plt.savefig("ba_gene_delta_star_boxplot.png", dpi=600)
plt.show()
