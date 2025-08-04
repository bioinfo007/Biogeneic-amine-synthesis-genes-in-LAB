#!/usr/bin/env python3
"""
plot_BA_gene_combinations.py
----------------------------

Description:
    Reads a gene presence/absence matrix and visualizes the distribution of
    biogenic amine (BA) gene co-occurrence combinations as a pie chart.

    Each genome's BA gene presence profile is converted into a combination label
    (e.g., "adc+tdc", "hdc", "No BA Genes"). The script counts all unique
    combinations across genomes, then plots their relative proportions.

Inputs:
    - gene_presence_absence_matrix.csv
        Must contain at least the following columns:
        adc, hdc, odc, tdc, agmatinase
        Values can be presence counts; nonzero values are treated as "present".

Outputs:
    - Pie chart of BA gene combinations (displayed interactively)
    - Can be easily modified to save as PNG/PDF

Usage:
    Place script in the same directory as gene_presence_absence_matrix.csv, then run:

        ./plot_BA_gene_combinations.py
    or:
        python3 plot_BA_gene_combinations.py

Requirements:
    - Python 3
    - pandas, numpy, matplotlib installed
      (e.g., pip install pandas numpy matplotlib)

Example:
    python3 plot_BA_gene_combinations.py

"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------
# Load data
# ---------------------------
df = pd.read_csv("gene_presence_absence_matrix.csv")

# Define BA genes
gene_cols = ["adc", "hdc", "odc", "tdc", "agmatinase"]

# Convert to binary (presence/absence)
df[gene_cols] = (df[gene_cols] > 0).astype(int)

# ---------------------------
# Create combination labels
# ---------------------------
def combination_label(row):
    genes = [gene for gene in gene_cols if row[gene] == 1]
    return "+".join(genes) if genes else "No BA Genes"

df["Combination"] = df.apply(combination_label, axis=1)

# Count combinations
combo_counts = df["Combination"].value_counts().reset_index()
combo_counts.columns = ["Combination", "Count"]

# ---------------------------
# Plot pie chart
# ---------------------------
colors = plt.cm.tab20(np.linspace(0, 1, len(combo_counts)))

def autopct_threshold(pct):
    return f"{pct:.1f}%" if pct >= 0 else ""

fig, ax = plt.subplots(figsize=(8, 8))
wedges, _, autotexts = ax.pie(
    combo_counts["Count"],
    labels=None,  # No direct labels on slices
    autopct=autopct_threshold,
    colors=colors,
    startangle=90
)

# Legend
ax.legend(
    wedges,
    combo_counts["Combination"],
    title="BA Gene Combination",
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    fontsize=12,
    title_fontsize=14
)

ax.set_title("BA Gene Co-occurrence Combinations Across Genomes", fontsize=16, fontweight="bold")

plt.tight_layout()
plt.show()
