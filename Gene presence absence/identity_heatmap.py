#!/usr/bin/env python3
"""
identity_heatmap.py
-------------------

Description:
    Plots an all-vs-all % identity heatmap with custom segmented gradient steps (0–30, 30–50, >50)
    using a provided identity matrix CSV.

Input:
    - identity_matrix.csv (rows and columns are sequence identifiers; values are percent identity 0–100)

Usage:
    ./identity_heatmap.py --input path/to/identity_matrix.csv [--output heatmap.svg]

Arguments:
    --input,  -i    Path to the identity matrix CSV. Required.  
                   First column is treated as index if present.
    --output, -o    Output filename for saved figure (SVG or PNG). Default: heatmap.svg

Requirements:
    - Python 3
    - pandas, seaborn, matplotlib, numpy installed
        e.g., pip install pandas seaborn matplotlib numpy

Example:
    python3 identity_heatmap.py -i identity_matrix.csv -o identity_heatmap.svg

Notes:
    - The script builds a custom colormap with 5% steps and highlights key breakpoints (0%, 30%, 50%, 100%).
    - Ensure the input file path is correct; relative or absolute paths are accepted.
"""

import argparse
import os
import sys

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm

def interpolate_colors(start_hex, end_hex, steps):
    import matplotlib.colors as mcolors
    start_rgb = np.array(mcolors.hex2color(start_hex))
    end_rgb = np.array(mcolors.hex2color(end_hex))
    return [mcolors.to_hex(start_rgb + (end_rgb - start_rgb) * i / (steps - 1)) for i in range(steps)]

def main():
    parser = argparse.ArgumentParser(description="Plot % identity heatmap with custom gradient.")
    parser.add_argument("-i", "--input", required=True, help="Identity matrix CSV file.")
    parser.add_argument("-o", "--output", default="heatmap.svg", help="Output figure filename (SVG/PNG).")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        sys.exit(f"ERROR: Input file not found: {args.input}")

    # Load matrix
    df = pd.read_csv(args.input, index_col=0)

    # Define color anchors per block
    block_colors = {
        '<30': ['#add8e6', '#0d47a1'],        # light blue → deep blue
        '30-50': ['#f4c29b', '#b2733f'],      # light tan → dark tan
        '>50': ['#d73027', '#7f0000'],         # light red → dark red
    }

    # Build segmented colormap with 5% steps
    steps1 = 7   # 0,5,...,30
    colors1 = interpolate_colors(block_colors['<30'][0], block_colors['<30'][1], steps1)
    positions1 = np.linspace(0, 0.3, steps1)

    steps2 = 5   # 30,35,...,50
    colors2 = interpolate_colors(block_colors['30-50'][0], block_colors['30-50'][1], steps2)
    positions2 = np.linspace(0.3, 0.5, steps2)

    steps3 = 11  # 50,55,...,100
    colors3 = interpolate_colors(block_colors['>50'][0], block_colors['>50'][1], steps3)
    positions3 = np.linspace(0.5, 1.0, steps3)

    # Combine avoiding duplicate boundaries
    positions = np.concatenate([positions1, positions2[1:], positions3[1:]])
    colors = colors1 + colors2[1:] + colors3[1:]

    custom_cmap = LinearSegmentedColormap.from_list("custom_5pct_cmap", list(zip(positions, colors)))

    norm = BoundaryNorm(boundaries=np.linspace(0, 100, 101), ncolors=custom_cmap.N)

    # Plot heatmap
    plt.figure(figsize=(15, 12))
    ax = sns.heatmap(
        df,
        cmap=custom_cmap,
        norm=norm,
        cbar_kws={'ticks': [0, 30, 50, 100], 'label': '% Identity'},
        square=True,
        linewidths=0.2,
        linecolor='white'
    )

    cbar = ax.collections[0].colorbar
    cbar.set_ticks([0, 30, 50, 100])
    cbar.set_ticklabels(['0%', '30%', '50%', '100%'])
    cbar.ax.tick_params(labelsize=12)

    plt.title("All-vs-All % Identity Heatmap with 5% Gradient Steps", fontsize=16)
    plt.xlabel("Sequences", fontsize=14, fontweight='bold', fontname='Arial')
    plt.ylabel("Sequences", fontsize=14, fontweight='bold', fontname='Arial')
    plt.xticks(rotation=90, fontsize=12, fontweight='bold', fontname='Arial')
    plt.yticks(rotation=0, fontsize=12, fontweight='bold', fontname='Arial')
    plt.tight_layout()

    # Save
    out = args.output
    ext = os.path.splitext(out)[1].lower()
    if ext not in ['.svg', '.png', '.pdf']:
        out = out + '.svg'
    plt.savefig(out, dpi=300, bbox_inches='tight')
    print(f"Saved heatmap to {out}")
    plt.show()

if __name__ == "__main__":
    main()
