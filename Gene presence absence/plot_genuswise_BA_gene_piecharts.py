#!/usr/bin/env python3
"""
plot_combination_by_genus.py
---------------------------

Description:
    For each genus, plots the co-occurrence combination of BA genes as a pie chart.
    Uses a presence/absence matrix with an "Organism" column and binary columns for
    adc, hdc, odc, tdc, agmatinase.

Inputs:
    - gene_presence_absence_matrix.csv
        Required columns: Organism, adc, hdc, odc, tdc, agmatinase

Outputs:
    - A multi-panel pie chart (one panel per genus) showing BA gene combinations, with
      sample size annotated in titles.
    - Displayed interactively and saved as:
        * BA_gene_combination_by_genus.pdf
        * BA_gene_combination_by_genus.png

Usage:
    Place the CSV in the same directory and run:
        ./plot_combination_by_genus.py
    or:
        python3 plot_combination_by_genus.py

Requirements:
    - Python 3
    - pandas, numpy, matplotlib installed (e.g., pip install pandas numpy matplotlib)
    Author:
    Aqib Javaid
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

# ---------------------------
# Load data
# ---------------------------
df = pd.read_csv("gene_presence_absence_matrix.csv", sep=",")  # adjust if needed

# Extract genus from 'Organism'
df["Genus"] = df["Organism"].str.split().str[0]

# Define gene columns (with Agmatinase)
gene_cols = ["adc", "hdc", "odc", "tdc", "agmatinase"]

# Binarize (1 if present, 0 if absent)
df[gene_cols] = (df[gene_cols] > 0).astype(int)

# Create combination label per genome
def combination_label(row):
    genes = [gene for gene in gene_cols if row[gene] == 1]
    return "+".join(genes) if genes else "No BA Genes"

df["Combination"] = df.apply(combination_label, axis=1)

# Get unique genera
unique_genera = sorted(df["Genus"].unique())
num_genera = len(unique_genera)

# Grid size for subplots
cols = 3
rows = math.ceil(num_genera / cols)

fig, axes = plt.subplots(rows, cols, figsize=(cols * 6, rows * 6))
axes = axes.flatten()

# Plot per genus
for i, genus in enumerate(unique_genera):
    ax = axes[i]
    sub = df[df["Genus"] == genus]
    counts = sub["Combination"].value_counts().reset_index()
    counts.columns = ["Combination", "Count"]
    total = counts["Count"].sum()

    colors = plt.cm.tab20(np.linspace(0, 1, len(counts)))

    wedges, _, autotexts = ax.pie(
        counts["Count"],
        labels=None,
        autopct="%1.1f%%",
        colors=colors,
        startangle=90,
        textprops={'fontsize': 10}
    )

    ax.legend(
        wedges,
        counts["Combination"],
        title="BA Genes",
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        fontsize=8,
        title_fontsize=10
    )

    ax.set_title(f"{genus} (n={total})", fontsize=14, fontweight="bold")

# Remove unused subplots
for j in range(len(unique_genera), len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()

# Save output
fig.savefig("BA_gene_combination_by_genus.pdf", dpi=300, bbox_inches='tight')
fig.savefig("BA_gene_combination_by_genus.png", dpi=300, bbox_inches='tight')

plt.show()
