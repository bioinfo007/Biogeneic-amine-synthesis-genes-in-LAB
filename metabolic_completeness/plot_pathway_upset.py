#!/usr/bin/env python3
"""
plot_pathway_upset.py
---------------------

Description:
    Reads pathway completeness matrix and plots an UpSet plot of pathway combinations
    (only genomes with at least one complete pathway are included).

Input:
    - pathway_completeness_matrix.csv
        Expected columns: Tyramine, Histamine, Putrescine_1, Putrescine_2, Putrescine_3, Agmatine
        Values should be binary (0/1) indicating completeness.

Output:
    - UpSet plot displayed and saved as:
        * pathway_combinations_upset.pdf
        * pathway_combinations_upset.png

Usage:
    Place this script alongside `pathway_completeness_matrix.csv` and run:

        ./plot_pathway_upset.py
    or:
        python3 plot_pathway_upset.py

Requirements:
    - Python 3
    - pandas, matplotlib, upsetplot installed
      (e.g., pip install pandas matplotlib upsetplot)
"""

import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_indicators
import sys
import os

# Input file
INPUT_FILE = "pathway_completeness_matrix.csv"
OUTPUT_PDF = "pathway_combinations_upset.pdf"
OUTPUT_PNG = "pathway_combinations_upset.png"

# Pathways to consider
PATHWAY_COLS = ["Tyramine", "Histamine", "Putrescine_1", "Putrescine_2", "Putrescine_3", "Agmatine"]

# --- Load data ---
if not os.path.isfile(INPUT_FILE):
    sys.exit(f"ERROR: Input file not found: {INPUT_FILE}")

df = pd.read_csv(INPUT_FILE, index_col=0)

# Check required columns
missing = [c for c in PATHWAY_COLS if c not in df.columns]
if missing:
    sys.exit(f"ERROR: Missing expected pathway columns in input: {missing}")

# Filter genomes with at least one complete pathway
df_filtered = df[df[PATHWAY_COLS].sum(axis=1) > 0]
if df_filtered.empty:
    sys.exit("ERROR: No genomes have any complete pathway; nothing to plot.")

# Boolean matrix for UpSet
bool_df = df_filtered[PATHWAY_COLS].astype(bool)

# Build UpSet data
data = from_indicators(bool_df, data=df_filtered)

# Plot
plt.figure(figsize=(10, 6))
upset = UpSet(
    data,
    show_counts=True,
    sort_by='cardinality',
    element_size=None,  # let it auto-scale
)

axes = upset.plot()

# Color customization: main intersection bars
if "intersections" in axes:
    for bar in axes["intersections"].patches:
        bar.set_facecolor("steelblue")
if "totals" in axes:
    for bar in axes["totals"].patches:
        bar.set_facecolor("steelblue")

plt.suptitle("Combination of Complete Biogenic Amine Pathways", fontsize=16, fontweight="bold")
plt.tight_layout(rect=[0, 0, 1, 0.95])

# Save outputs
plt.savefig(OUTPUT_PDF, dpi=300, bbox_inches='tight')
plt.savefig(OUTPUT_PNG, dpi=300, bbox_inches='tight')

plt.show()
