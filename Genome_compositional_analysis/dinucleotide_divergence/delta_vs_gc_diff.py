"""
Scatter plot of GC% difference (gene vs genome) versus Karlin δ* for biogenic amine (BA) genes.
Colored by gene type with threshold lines to highlight compositional bias cutoffs.

Input:
- 'combine_compositional_result.csv' with columns: 'gc_diff', 'Delta_Star', 'Gene_Name'

Output:
- Saves figure as 'delta_vs_gc_diff.png'

Usage:
- Make sure the input CSV is in the working directory.
- Run this script to generate and save the scatter plot.
Author:
    Aqib Javaid
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("combine_compositional_result.csv")

# Plot GC diff vs δ*
plt.figure(figsize=(8, 6))
sns.set(style="whitegrid")

scatter = sns.scatterplot(
    data=df,
    x="gc_diff",
    y="Delta_Star",
    hue="Gene_Name",
    palette="Set1",
    s=80,
    alpha=0.8,
    edgecolor="black"
)

# Add horizontal (δ*) threshold line (adjust threshold, if needed)
plt.axhline(y=0.145, color='gray', linestyle='--', linewidth=1.5, label=r"δ* threshold (0.145)")

# Add vertical (GC%) threshold line (adjust threshold, if needed)
plt.axvline(x=3.0, color='black', linestyle='--', linewidth=1.5, label="GC% diff threshold (3%)")

# Customize labels and title
plt.xlabel("GC% Difference (Gene - Genome)", fontsize=12)
plt.ylabel(r"δ* (Karlin Dinucleotide Signature Difference)", fontsize=12)
plt.title("δ* vs. GC% Difference in BA Genes", fontsize=14, weight='bold')

# Add legend outside plot
plt.legend(title="Gene", fontsize=10, title_fontsize=11, loc='upper left', bbox_to_anchor=(1, 1))

plt.tight_layout()
plt.savefig("delta_vs_gc_diff.png", dpi=600, bbox_inches='tight')
plt.show()
