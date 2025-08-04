#!/usr/bin/env python3
"""
plot_BA_gene_presence.py
------------------------

Description:
    Computes and plots the percentage of genomes per genus carrying each biogenic amine (BA)
    gene from a presence/absence matrix. Saves the per-genus presence percentages and
    generates publication-ready bar plots.

Inputs:
    - gene_presence_absence_matrix.csv : CSV with columns including "Organism", and binary
      presence/absence columns for genes: adc, hdc, odc, tdc, agmatinase.
        * Genome IDs in "Organism" are expected to start with genus/species name, e.g., "Lactobacillus xyz".
        * Presence columns should be 0/1.

Outputs:
    - BA_gene_presence_percentage_per_genus.csv : percentage presence of each gene per genus.
    - BA_gene_presence_per_genus.pdf / .png : bar plot of gene presence percentages by genus.

Usage:
    Place this script in the directory containing `gene_presence_absence_matrix.csv`, then run:

        ./plot_BA_gene_presence.py

    Or explicitly:

        python3 plot_BA_gene_presence.py

Requirements:
    - Python 3
    - pandas, numpy, matplotlib installed (e.g., via `pip install pandas numpy matplotlib`)
    - Font "Arial" available (optional; falls back if not present)

Example:
    python3 plot_BA_gene_presence.py
    Author:
    Aqib Javaid

"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------
# Font / style settings
# ---------------------------
font_family = 'Arial'
title_fontsize = 18
title_fontweight = 'bold'
label_fontsize = 16
label_fontweight = 'bold'
tick_fontsize = 14
tick_fontweight = 'bold'
legend_fontsize = 14
legend_title_fontsize = 16
annotation_fontsize = 12

plt.rcParams.update({
    "font.family": font_family,
})

# ---------------------------
# Load data
# ---------------------------
df = pd.read_csv("gene_presence_absence_matrix.csv")  # expects comma-separated

# Extract Genus from organism name column
df["Genus"] = df["Organism"].str.split().str[0]

# Define genes of interest
gene_cols = ["adc", "hdc", "odc", "tdc", "agmatinase"]

# Ensure presence/absence are integers 0/1
df[gene_cols] = df[gene_cols].apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)

# Compute percentage presence per genus
gene_presence = df.groupby("Genus")[gene_cols].mean() * 100  # percent

# Save the table
gene_presence.reset_index().to_csv("BA_gene_presence_percentage_per_genus.csv", index=False)

# Get total genomes per genus for annotation
total_genomes = df.groupby("Genus").size()

# ---------------------------
# Plotting
# ---------------------------
genera = gene_presence.index.tolist()
num_genera = len(genera)
num_genes = len(gene_cols)

# Bar positioning
bar_width = 0.15
x = np.arange(num_genera)
offsets = (np.arange(num_genes) - (num_genes - 1)/2) * bar_width  # centered

# Colors for genes
presence_colors = ["#4ee90a", "#0acbf3", "#e71f10", "#0000ff", "#6e03f8"]

fig, ax = plt.subplots(figsize=(16, 7))

for i, gene in enumerate(gene_cols):
    ax.bar(
        x + offsets[i],
        gene_presence[gene],
        width=bar_width,
        label=gene.replace("_", " ").title(),
        color=presence_colors[i],
        edgecolor='black',
        linewidth=0.5,
    )

# X-axis
ax.set_xticks(x)
ax.set_xticklabels(
    genera,
    rotation=45,
    ha='right',
    fontsize=tick_fontsize,
    fontweight=tick_fontweight,
)

# Labels and title
ax.set_ylabel("Percentage of Genomes with Gene (%)", fontsize=label_fontsize, fontweight=label_fontweight)
ax.set_xlabel("Genus", fontsize=label_fontsize, fontweight=label_fontweight)
ax.set_title("Proportion of Genomes with BA Genes per Genus", fontsize=title_fontsize, fontweight=title_fontweight)

# Legend
legend = ax.legend(title="BA Genes", fontsize=legend_fontsize, title_fontsize=legend_title_fontsize,
                   loc='upper left', bbox_to_anchor=(1.02, 1))
ax.add_artist(legend)

# Annotate sample sizes below x-axis
for idx, genus in enumerate(genera):
    ax.text(
        idx,
        -5,  # a bit below zero
        f"n={total_genomes[genus]}",
        ha='center',
        va='top',
        rotation=45,
        fontsize=annotation_fontsize,
        clip_on=False
    )

# Adjust limits and layout
ax.set_ylim(0, 105)
plt.tight_layout()

# Save figure
fig.savefig("BA_gene_presence_per_genus.pdf", dpi=300, bbox_inches='tight')
fig.savefig("BA_gene_presence_per_genus.png", dpi=300, bbox_inches='tight')

plt.show()
