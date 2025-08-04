#!/usr/bin/env python3
"""
compute_gc_skew.py
--------------------

Description:
    For each GenBank/GBFF file in a folder, finds CDS features whose locus_tag
    matches any in a provided list, and compares GC content and GC skew of the
    gene vs the whole genome.

Input:
    - Folder of GenBank (.gbk/.gbff) files (e.g., Prokka outputs)
    - Plain text file of BA gene locus_tag IDs (one per line)

Output:
    - TSV with per-gene and per-genome GC% and GC skew, plus differences.

Usage:
    ./gc_skew_locus_tag.py -g path/to/gbk_folder -l locus_tags.txt -o output.tsv

Requirements:
    - Python 3
    - Biopython
    - pandas

Example:
    ./compute_gc_skew.py -g prokka_gbks/ -l adc_locus_tags.txt -o adc_gc_comparison.tsv
    Author:
    Aqib Javaid
"""
import os
from Bio import SeqIO
import argparse
import pandas as pd

def gc_content(seq):
    """Calculate GC content (%) and GC skew of a sequence."""
    seq = seq.upper()
    g = seq.count('G')
    c = seq.count('C')
    total = len(seq)
    gc_pct = ((g + c) / total) * 100 if total > 0 else 0
    gc_skew = (g - c) / (g + c) if (g + c) > 0 else 0
    return gc_pct, gc_skew

def extract_gene_by_locus_tag(gbk_file, target_tags):
    """Find a gene with matching locus_tag in the GBK and calculate GC metrics."""
    for record in SeqIO.parse(gbk_file, "genbank"):
        genome_seq = record.seq
        genome_gc, genome_skew = gc_content(str(genome_seq))
        for feature in record.features:
            if feature.type == "CDS" and "locus_tag" in feature.qualifiers:
                tag = feature.qualifiers["locus_tag"][0]
                if tag in target_tags:
                    gene_seq = feature.extract(genome_seq)
                    gene_gc, gene_skew = gc_content(str(gene_seq))
                    return tag, gene_gc, gene_skew, genome_gc, genome_skew
    return None, None, None, None, None

def main(gbk_folder, locus_tag_file, output_file):
    # Load locus_tag list
    with open(locus_tag_file) as f:
        locus_tags = set(line.strip() for line in f if line.strip())

    results = []

    for fname in os.listdir(gbk_folder):
        if fname.endswith(".gbk") or fname.endswith(".gbff"):
            gbk_path = os.path.join(gbk_folder, fname)
            genome_name = os.path.splitext(fname)[0]
            tag, gene_gc, gene_skew, genome_gc, genome_skew = extract_gene_by_locus_tag(gbk_path, locus_tags)

            if tag:
                results.append({
                    "genome": genome_name,
                    "locus_tag": tag,
                    "gene_gc_percent": round(gene_gc, 3),
                    "genome_gc_percent": round(genome_gc, 3),
                    "gene_gc_skew": round(gene_skew, 4),
                    "genome_gc_skew": round(genome_skew, 4),
                    "gc_diff": round(gene_gc - genome_gc, 3),
                    "skew_diff": round(gene_skew - genome_skew, 4)
                })

    df = pd.DataFrame(results)
    df.to_csv(output_file, sep='\t', index=False)

    print(f"\n✅ Results saved to: {output_file}")
    print(f"Analyzed {len(results)} genomes with matching locus_tags.\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare GC content and GC skew for genes using locus_tag IDs.")
    parser.add_argument("-g", "--gbk_folder", required=True, help="Folder with Prokka-annotated .gbk files")
    parser.add_argument("-l", "--locus_tags", required=True, help="Text file with locus_tag list (one per line)")
    parser.add_argument("-o", "--output", default="gc_skew_locus_tag_output.tsv", help="Output TSV file name")
    args = parser.parse_args()

    main(args.gbk_folder, args.locus_tags, args.output)
