#!/usr/bin/env python3
"""
build_presence_absence.py
-------------------------

Description:
    Builds a gene presence/absence matrix from a genome list and per-gene FASTA files.

Usage:
    ./build_presence_absence.py --genome-list genomes.txt --fasta-dir extracted_BA_genes --output BA_gene_matrix.csv

Arguments:
    --genome-list, -g   Path to genome list file (one genome ID per line, e.g., SAMN037031)
    --fasta-dir,  -f    Path to directory containing gene FASTA files
    --output,     -o    Output CSV filename (default: gene_presence_absence_matrix.csv)

Requirements:
    - Python 3
    - Biopython installed (`pip install biopython`)
    - Genome IDs in genomes.txt must match the prefix before the first underscore in FASTA headers
    - FASTA files should be named as per the mapping in the `genes` dictionary inside the script

Example:
    ./build_presence_absence.py \
        --genome-list genomes.txt \
        --fasta-dir extracted_BA_genes \
        --output BA_gene_matrix.csv

Author:
    Aqib Javaid
"""

import argparse
import csv
import os
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(
        description="Build a gene presence/absence matrix from genome list and gene FASTA files."
    )
    parser.add_argument(
        "--genome-list", "-g", required=True,
        help="Path to genome list file (one genome ID per line)"
    )
    parser.add_argument(
        "--fasta-dir", "-f", required=True,
        help="Path to directory containing FASTA files (one per gene)"
    )
    parser.add_argument(
        "--output", "-o", default="gene_presence_absence_matrix.csv",
        help="Output CSV file name"
    )

    args = parser.parse_args()

    # Define your gene-to-file mapping here
    genes = {
        "adc": "adc.faa.pep",
        "hdc": "hdc.faa.pep",
        "agmatinase": "agmatinase.faa.pep",
        "odc": "odc.pep",
        "tdc": "tdc.pep"
    }

    # Load genome list
    if not os.path.isfile(args.genome_list):
        parser.error(f"Genome list file not found: {args.genome_list}")
    with open(args.genome_list) as f:
        all_genomes = sorted(line.strip() for line in f if line.strip())

    # Initialize presence matrix
    presence = {genome: {gene: 0 for gene in genes} for genome in all_genomes}

    # Scan each gene FASTA file
    for gene, fasta_name in genes.items():
        fasta_path = os.path.join(args.fasta_dir, fasta_name)
        if not os.path.isfile(fasta_path):
            print(f"WARNING: FASTA file not found for {gene}: {fasta_path}")
            continue
        for record in SeqIO.parse(fasta_path, "fasta"):
            samn_id = record.id.split("_")[0]  # Prefix before first underscore
            if samn_id in presence:
                presence[samn_id][gene] = 1

    # Write output CSV
    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Genome"] + list(genes.keys()))
        for genome in all_genomes:
            row = [genome] + [presence[genome][gene] for gene in genes]
            writer.writerow(row)

    print(f"Matrix saved to: {args.output}")

if __name__ == "__main__":
    main()
