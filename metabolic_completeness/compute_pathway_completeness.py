#!/usr/bin/env python3
"""
compute_pathway_completeness.py
-------------------------------

Description:
    Given a gene presence/absence matrix, evaluates completeness of predefined
    metabolic pathways per genome. Each pathway is defined as a set of required
    genes; a genome is marked complete (1) if all required genes are present.

Input:
    - gene presence/absence CSV (rows = genomes, columns = genes; values 0/1)

Output:
    - pathway_completeness_matrix.csv: binary matrix of pathway completeness per genome

Usage:
    ./compute_pathway_completeness.py --input gene_presence_absence_matrix.csv \
        --output pathway_completeness_matrix.csv

    Optional: supply a custom pathway definition file (JSON) with --pathways

Defaults:
    Input: gene_presence_absence_matrix.csv
    Output: pathway_completeness_matrix.csv
    Pathways: built-in definitions (Tyramine, Histamine, Putrescine variants, Agmatine)

Requirements:
    - Python 3
    - pandas installed (`pip install pandas`)
"""

import argparse
import json
import pandas as pd
import sys
import os

# Built-in pathway definitions
DEFAULT_PATHWAYS = {
    'Tyramine': ['tdc'],
    'Histamine': ['hdc'],
    'Putrescine_1': ['adi', 'otc', 'odc'],      # ADI + OTC + ODC
    'Putrescine_2': ['adc', 'agdi', 'ptc'],     # ADC + AgDI + PTC
    'Putrescine_3': ['adc', 'agmatinase'],      # ADC + Agmatinase
    'Agmatine': ['adc']
}

def load_pathways(pathway_file):
    with open(pathway_file) as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError("Pathway file must be a JSON object mapping pathway names to gene lists.")
    return data

def main():
    parser = argparse.ArgumentParser(description="Compute pathway completeness from gene presence/absence matrix.")
    parser.add_argument(
        "--input", "-i", default="gene_presence_absence_matrix.csv",
        help="Input gene presence/absence CSV (genomes x genes, values 0/1)."
    )
    parser.add_argument(
        "--output", "-o", default="pathway_completeness_matrix.csv",
        help="Output CSV for pathway completeness."
    )
    parser.add_argument(
        "--pathways", "-p", default=None,
        help="Optional JSON file defining pathways. Format: {\"PathwayName\": [\"gene1\",\"gene2\", ...], ...}"
    )
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        sys.exit(f"ERROR: Input file not found: {args.input}")

    # Load gene matrix
    df_genes = pd.read_csv(args.input, index_col=0)
    # Ensure values are binary 0/1
    df_genes = (df_genes > 0).astype(int)

    # Load pathways
    if args.pathways:
        if not os.path.isfile(args.pathways):
            sys.exit(f"ERROR: Pathway definition file not found: {args.pathways}")
        try:
            pathways = load_pathways(args.pathways)
        except Exception as e:
            sys.exit(f"ERROR loading pathways JSON: {e}")
    else:
        pathways = DEFAULT_PATHWAYS

    df_pathways = pd.DataFrame(index=df_genes.index)

    for pathway_name, gene_list in pathways.items():
        missing = [g for g in gene_list if g not in df_genes.columns]
        if missing:
            print(f"Warning: genes {missing} not found in matrix for pathway '{pathway_name}'. Assuming absence.", file=sys.stderr)
        genes_present = [g for g in gene_list if g in df_genes.columns]
        if genes_present:
            completeness = df_genes[genes_present].all(axis=1).astype(int)
        else:
            completeness = pd.Series(0, index=df_genes.index)
        df_pathways[pathway_name] = completeness

    df_pathways.to_csv(args.output)
    print(f"Pathway completeness matrix saved to '{args.output}'")

if __name__ == "__main__":
    main()
