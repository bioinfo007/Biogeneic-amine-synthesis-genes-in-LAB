#!/usr/bin/env python3
"""
gc_skew_plotter.py
------------------

Description:
    For each gene-specific GC/skew TSV (expected format from prior GC analysis),
    generates:
      A. Violin+swarm plots comparing gene vs genome GC content and skew, with difference histograms.
      B. 2D scatter of GC content vs skew for gene and genome.

Expected input files:
    <input_dir>/<GENE>_GC.tsv
    Each should contain at least the columns:
        genome, gene_gc_percent, genome_gc_percent, gene_gc_skew, genome_gc_skew

Defaults:
    input_dir: ./gc_data
    output_dir: ./plots
    genes: ADC, ODC, TDC, HDC, agmatinase

Usage:
    ./gc_skew_plotter.py \
        --input-dir path/to/gc_data \
        --output-dir path/to/plots \
        --genes ADC ODC TDC HDC agmatinase

Requirements:
    - Python 3
    - pandas, seaborn, matplotlib, numpy

Example:
    python3 gc_skew_plotter.py --input-dir ./gc_data --output-dir ./plots
    Author:
    Aqib Javaid
"""

import os
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Plot gene vs genome GC content and skew.")
    p.add_argument(
        "--input-dir", "-i", default="./gc_data",
        help="Directory containing per-gene *_GC.tsv files."
    )
    p.add_argument(
        "--output-dir", "-o", default="./plots",
        help="Directory to save generated plots."
    )
    p.add_argument(
        "--genes", "-g", nargs="+", default=["ADC", "ODC", "TDC", "HDC", "agmatinase"],
        help="List of gene prefixes (will look for <gene>_GC.tsv)."
    )
    return p.parse_args()

def plot_gene_gc_skew(df, gene_name, output_dir, params):
    # Prepare data for GC content violin + box
    gc_melt = pd.melt(df, id_vars=["genome"], value_vars=["gene_gc_percent", "genome_gc_percent"],
                      var_name="Region", value_name="GC_Content")
    gc_melt["Region"] = gc_melt["Region"].map({"gene_gc_percent": "Gene", "genome_gc_percent": "Genome"})

    # Prepare data for GC skew violin + box
    skew_melt = pd.melt(df, id_vars=["genome"], value_vars=["gene_gc_skew", "genome_gc_skew"],
                        var_name="Region", value_name="GC_Skew")
    skew_melt["Region"] = skew_melt["Region"].map({"gene_gc_skew": "Gene", "genome_gc_skew": "Genome"})

    # Delta dataframe for difference plots
    df["gc_diff"] = df["gene_gc_percent"] - df["genome_gc_percent"]
    df["skew_diff"] = df["gene_gc_skew"] - df["genome_gc_skew"]

    # Start plotting
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel A: GC content violin + box + swarm
    sns.violinplot(x="Region", y="GC_Content", data=gc_melt, palette=params["palette"], inner="box", ax=axes[0,0])
    sns.swarmplot(x="Region", y="GC_Content", data=gc_melt, color=params["swarm_color"], size=params["swarm_size"], ax=axes[0,0])
    axes[0,0].set_title(f"{gene_name} - GC Content", fontsize=params["fontsize_title"], fontfamily=params["fontfamily"])
    axes[0,0].set_ylabel("GC Content (%)", fontsize=params["fontsize_axes_labels"], fontfamily=params["fontfamily"])
    axes[0,0].tick_params(axis='both', labelsize=params["fontsize_ticks"])

    # Panel B: GC skew violin + box + swarm with zero line
    sns.violinplot(x="Region", y="GC_Skew", data=skew_melt, palette=params["palette"], inner="box", ax=axes[0,1])
    sns.swarmplot(x="Region", y="GC_Skew", data=skew_melt, color=params["swarm_color"], size=params["swarm_size"], ax=axes[0,1])
    axes[0,1].axhline(0, ls='--', c='grey')
    axes[0,1].set_title(f"{gene_name} - GC Skew", fontsize=params["fontsize_title"], fontfamily=params["fontfamily"])
    axes[0,1].set_ylabel("GC Skew", fontsize=params["fontsize_axes_labels"], fontfamily=params["fontfamily"])
    axes[0,1].tick_params(axis='both', labelsize=params["fontsize_ticks"])

    # Panel C: Histogram of GC diff
    sns.histplot(df["gc_diff"], bins=params["hist_bins"], color=params["palette"]["Gene"], alpha=params["hist_alpha"], label="GC Content Diff", ax=axes[1,0])
    axes[1,0].axvline(0, ls='--', c='grey')
    axes[1,0].set_title(f"{gene_name} - GC Content Difference (Gene - Genome)", fontsize=params["fontsize_title"], fontfamily=params["fontfamily"])
    axes[1,0].set_xlabel("GC Content Difference (%)", fontsize=params["fontsize_axes_labels"], fontfamily=params["fontfamily"])
    axes[1,0].legend()
    axes[1,0].tick_params(axis='both', labelsize=params["fontsize_ticks"])

    # Panel D: Histogram of Skew diff
    sns.histplot(df["skew_diff"], bins=params["hist_bins"], color=params["palette"]["Genome"], alpha=params["hist_alpha"], label="GC Skew Diff", ax=axes[1,1])
    axes[1,1].axvline(0, ls='--', c='grey')
    axes[1,1].set_title(f"{gene_name} - GC Skew Difference (Gene - Genome)", fontsize=params["fontsize_title"], fontfamily=params["fontfamily"])
    axes[1,1].set_xlabel("GC Skew Difference", fontsize=params["fontsize_axes_labels"], fontfamily=params["fontfamily"])
    axes[1,1].legend()
    axes[1,1].tick_params(axis='both', labelsize=params["fontsize_ticks"])

    plt.suptitle(f"GC Content and Skew Analysis for {gene_name}", fontsize=params["fontsize_title"] + 2, fontfamily=params["fontfamily"])
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out_path = os.path.join(params["output_dir"], f"{gene_name}_GC_Skew_analysis.svg")
    plt.savefig(out_path)
    plt.close()

def plot_2d_gc_skew(df, gene_name, params):
    plt.figure(figsize=(8,6))
    plt.scatter(df["genome_gc_percent"], df["genome_gc_skew"], label="Genome", s=params["scatter_size"], alpha=params["scatter_alpha"], c=params["palette"]["Genome"])
    plt.scatter(df["gene_gc_percent"], df["gene_gc_skew"], label="Gene", s=params["scatter_size"], alpha=params["scatter_alpha"], c=params["palette"]["Gene"])
    plt.xlabel("GC Content (%)", fontsize=params["fontsize_axes_labels"], fontfamily=params["fontfamily"])
    plt.ylabel("GC Skew", fontsize=params["fontsize_axes_labels"], fontfamily=params["fontfamily"])
    plt.title(f"GC Content vs GC Skew: {gene_name} Gene vs Genome", fontsize=params["fontsize_title"], fontfamily=params["fontfamily"])
    plt.legend()
    if params["grid"]:
        plt.grid(True)
    plt.tight_layout()
    out_path = os.path.join(params["output_dir"], f"{gene_name}_GCvsSkew_scatter.svg")
    plt.savefig(out_path)
    plt.close()

def main():
    args = parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir
    genes = args.genes

    os.makedirs(output_dir, exist_ok=True)

    # Plotting style and params
    sns.set(style="whitegrid")
    plot_params = {
        "palette": {"Gene": "#52aaec", "Genome": "#7af27a"},
        "fontsize_title": 16,
        "fontsize_axes_labels": 14,
        "fontsize_ticks": 12,
        "fontfamily": "Arial",
        "text_color": "black",
        "fontweight": "bold",
        "swarm_color": ".3",
        "swarm_size": 5,
        "hist_alpha": 0.7,
        "hist_bins": 25,
        "scatter_alpha": 0.7,
        "scatter_size": 40,
        "grid": True,
        "output_dir": output_dir
    }

    for gene in genes:
        # Expect filenames like ADC_GC.tsv or ADC.tsv? adapt:
        # trying both <gene>_GC.tsv and <gene>.tsv
        candidates = [f"{gene}_GC.tsv", f"{gene}.tsv"]
        file_path = None
        for cand in candidates:
            p = os.path.join(input_dir, cand)
            if os.path.exists(p):
                file_path = p
                break
        if not file_path:
            print(f"Warning: no file found for {gene} (looked for {candidates}), skipping.")
            continue

        df_gene = pd.read_csv(file_path, sep="\t")
        expected_cols = {"genome", "gene_gc_percent", "genome_gc_percent", "gene_gc_skew", "genome_gc_skew"}
        if not expected_cols.issubset(df_gene.columns):
            print(f"ERROR: {file_path} missing expected columns: {expected_cols - set(df_gene.columns)}. Skipping.")
            continue

        print(f"Processing {gene} ({df_gene.shape[0]} genomes)")

        plot_gene_gc_skew(df_gene, gene, output_dir, plot_params)
        plot_2d_gc_skew(df_gene, gene, plot_params)

    print("All plots saved in:", output_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GC content/skew plotting for genes vs genomes.")
    parser.add_argument("--input-dir", "-i", default="./gc_data", help="Directory with per-gene GC TSVs.")
    parser.add_argument("--output-dir", "-o", default="./plots", help="Directory to write plots.")
    parser.add_argument("--genes", "-g", nargs="+", default=["ADC", "ODC", "TDC", "HDC", "agmatinase"],
                        help="Gene prefixes to process (will look for <gene>_GC.tsv or <gene>.tsv).")
    args = parser.parse_args()

    # expose args to main
    def parse_args():
        return args

    main()
