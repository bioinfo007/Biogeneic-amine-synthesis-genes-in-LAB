#!/usr/bin/env python3
"""
plot_temporal_BA_genes.py
-------------------------

Description:
    Plots the temporal distribution (year-wise counts) of biogenic amine biosynthesis genes
    across LAB genomes. Each gene is shown in its own facet with:
      - Line plot of counts per year
      - First detection year marked
      - Major ticks every 20 years and minor ticks every year

Input:
    - gene_presence_absence_matrix.csv
        Required columns:
            * Genome
            * Date (containing a year like 1998, 2020, etc.)
            * Presence/absence columns for genes: adc, hdc, agmatinase, odc, tdc

Output:
    - SVG figure: BA_gene_individual_panels_with_major_minor_ticks.svg

Usage:
    Place this script in the same directory as the input CSV and run:
        ./plot_temporal_BA_genes.py
    or:
        python3 plot_temporal_BA_genes.py

Requirements:
    - Python 3
    - pandas, matplotlib, seaborn installed
        Install via: pip install pandas matplotlib seaborn

Notes:
    - Dates are scanned for a 4-digit year between 1900 and 2099.
    - Genes with zero counts across all years will still appear but may have empty trends.
    Author:
    Aqib Javaid
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import math
import sys
import os

# ---- CONFIG ----
INPUT_FILE = "gene_presence_absence_matrix.csv"
OUTPUT_SVG = "BA_gene_individual_panels_with_major_minor_ticks.svg"
GENES = ['adc', 'hdc', 'agmatinase', 'odc', 'tdc']
MIN_YEAR_PADDING = 2  # padding before/after observed years for xlim
MAJOR_TICK_INTERVAL = 20  # years between major ticks

# ---- Helper functions ----
def extract_year(date):
    if pd.isna(date):
        return None
    match = re.search(r'(19|20)\d{2}', str(date))
    return int(match.group()) if match else None

# ---- Load data ----
if not os.path.isfile(INPUT_FILE):
    sys.exit(f"ERROR: Input file not found: {INPUT_FILE}")

df = pd.read_csv(INPUT_FILE, sep=",", dtype=str)  # read all as str to safely parse

# Check required columns
required_base_cols = ["Genome", "Date"]
missing_req = [c for c in required_base_cols if c not in df.columns]
if missing_req:
    sys.exit(f"ERROR: Missing required columns in input: {missing_req}")

# Ensure gene columns exist; if missing, warn and drop
existing_genes = [g for g in GENES if g in df.columns]
missing_genes = [g for g in GENES if g not in df.columns]
if missing_genes:
    print(f"WARNING: The following genes are missing from input and will be skipped: {missing_genes}", file=sys.stderr)
if not existing_genes:
    sys.exit("ERROR: None of the expected gene columns present. Aborting.")

# Extract year
df['Year'] = df['Date'].apply(extract_year)
df = df.dropna(subset=['Year'])
df['Year'] = df['Year'].astype(int)

# Binarize gene presence for existing genes: treat any non-zero/non-empty as presence
for gene in existing_genes:
    # convert to numeric if possible, non-numeric treated as 0/absent
    df[gene] = pd.to_numeric(df[gene], errors='coerce').fillna(0).astype(int)
    df[gene] = (df[gene] > 0).astype(int)

# Reshape to long format for plotting
long_df = df.melt(id_vars=['Genome', 'Year'], value_vars=existing_genes,
                  var_name='Gene', value_name='Presence')

# Summarize counts per year per gene
summary = (
    long_df.groupby(['Year', 'Gene'], observed=True)
           .agg(gene_count=('Presence', 'sum'))
           .reset_index()
)

# First detection year per gene
first_appearance = (
    summary[summary['gene_count'] > 0]
    .groupby('Gene')['Year']
    .min()
    .to_dict()
)

# ---- Plotting ----
sns.set(style="whitegrid", context="talk", font_scale=1.0)
palette = {
    'adc': '#e41a1c', 'hdc': '#377eb8',
    'agmatinase': '#4daf4a', 'odc': '#984ea3', 'tdc': '#ff7f00'
}

# Prepare FacetGrid: ensure gene order matches existing_genes
g = sns.FacetGrid(summary, col="Gene", col_wrap=3, height=4, aspect=1.2, sharey=False,
                  col_order=existing_genes)

for ax, gene in zip(g.axes.flatten(), existing_genes):
    subset = summary[summary['Gene'] == gene]
    if subset.empty:
        ax.text(0.5, 0.5, f"No data for {gene}", ha='center', va='center')
        ax.set_title(gene.upper(), fontsize=13, weight='bold')
        ax.set_xlabel("Year")
        ax.set_ylabel("Gene Count")
        continue

    sns.lineplot(data=subset, x="Year", y="gene_count", color=palette.get(gene, 'gray'),
                 ax=ax, linewidth=2.2, marker="o")

    # Mark first appearance
    if gene in first_appearance:
        year = first_appearance[gene]
        ax.axvline(year, color='gray', linestyle='--', linewidth=1.5)
        ylim = ax.get_ylim()
        y_text = ylim[1] * 0.8 if ylim[1] > 0 else 0
        ax.text(year + 0.2, y_text, f"First: {year}", rotation=90,
                va='center', ha='left', fontsize=9, color='gray')

    # Axis limits with padding
    min_year = subset['Year'].min()
    max_year = subset['Year'].max()
    ax.set_xlim(min_year - MIN_YEAR_PADDING, max_year + MIN_YEAR_PADDING)

    # Major ticks every MAJOR_TICK_INTERVAL years
    start_tick = (min_year // MAJOR_TICK_INTERVAL) * MAJOR_TICK_INTERVAL
    end_tick = ((max_year // MAJOR_TICK_INTERVAL) + 1) * MAJOR_TICK_INTERVAL
    major_ticks = list(range(start_tick, end_tick + 1, MAJOR_TICK_INTERVAL))
    ax.set_xticks(major_ticks)

    # Minor ticks every year
    minor_ticks = list(range(start_tick, end_tick + 1, 1))
    ax.set_xticks(minor_ticks, minor=True)

    # Tick customization
    ax.tick_params(axis='x', which='major', length=7, width=1.5)
    ax.tick_params(axis='x', which='minor', length=3, width=0.8)
    ax.set_title(gene.upper(), fontsize=13, weight='bold')
    ax.set_xlabel("Year", fontsize=11)
    ax.set_ylabel("Gene Count", fontsize=11)
    for label in ax.get_xticklabels():
        label.set_rotation(45)
        label.set_ha('right')

# Clean up unused facets if any (FacetGrid handles this automatically)

plt.subplots_adjust(top=0.9)
g.fig.suptitle("Temporal Distribution of Biogenic Amine Biosynthesis Genes in LAB Genomes",
               fontsize=16, weight='bold')

# Save high-resolution SVG and PNG
g.savefig(OUTPUT_SVG, format='svg', dpi=600)
png_out = OUTPUT_SVG.replace('.svg', '.png')
g.savefig(png_out, dpi=300)

plt.show()
