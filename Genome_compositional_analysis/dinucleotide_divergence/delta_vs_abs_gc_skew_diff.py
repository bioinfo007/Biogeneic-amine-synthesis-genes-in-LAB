"""
Scatter plot of absolute GC skew difference (gene vs genome) versus Karlin δ* for biogenic amine (BA) genes.
Colored by gene type with cutoff lines indicating compositional bias thresholds.

Input:
- 'combine_compositional_result.csv' with columns: 'skew_diff', 'Delta_Star', 'Gene_Name'

Output:
- Saves figure as 'delta_vs_abs_gc_skew_diff.png'

Usage:
- Ensure input CSV is in the working directory.
- Run this script to generate and save the scatter plot.
Author:
    Aqib Javaid
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("combine_compositional_result.csv")

# Use absolute value of skew_diff for x-axis
df['abs_skew_diff'] = df['skew_diff'].abs()

# Thresholds
delta_cutoff = 0.145      # δ* cutoff for HGT
skew_cutoff = 0.15        # absolute GC skew difference cutoff

# Plot abs GC skew diff vs δ*
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=df,
    x="abs_skew_diff",
    y="Delta_Star",
    hue="Gene_Name",
    palette="Set1",
    s=80,
    alpha=0.8,
    edgecolor="black"
)

# Add cutoff lines
plt.axhline(delta_cutoff, linestyle='--', color='red', lw=2, label=f'δ* cutoff ({delta_cutoff})')
plt.axvline(skew_cutoff, linestyle='--', color='blue', lw=2, label=f'GC skew diff cutoff ({skew_cutoff})')

plt.xlabel("Absolute GC Skew Difference (|Gene - Genome|)", fontsize=12)
plt.ylabel(r"δ* (Karlin Dinucleotide Signature Difference)", fontsize=12)
plt.title("δ* vs. Absolute GC Skew Difference in BA Genes", fontsize=14, weight='bold')
plt.legend()
plt.tight_layout()
plt.savefig("delta_vs_abs_gc_skew_diff.png", dpi=600)
plt.show()
