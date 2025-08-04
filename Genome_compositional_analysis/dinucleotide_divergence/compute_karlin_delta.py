#!/usr/bin/env python3
"""
Compute Karlin δ* (delta star) metric comparing dinucleotide relative abundances
of genes versus their genome background sequences.

Input:
- Mapping CSV with columns: Genome, Gene_ID, (optional) Gene_Name
- Genome directory containing paired genome FASTA (.fna) and annotation GFF (.gff) files
  named consistently by genome ID.

For each gene:
- Extract gene CDS sequence from GFF and genome FASTA
- Extract longest contig sequence from genome FASTA as background
- Compute Karlin rho (relative dinucleotide abundance) for gene and genome background
- Calculate Karlin δ* as average absolute difference between gene and genome rho vectors
- Filter genes by minimum length (default 30 bp)
- Save results to CSV with Genome, Gene_ID, Gene_Name, Gene_Length, Delta_Star

Usage example:
$ python karlin_delta_star.py -m mapping.csv -g genome_dir -o output.csv

Requirements:
- Biopython
- pandas
- argparse
- logging

Note:
- GFF file must have CDS features with IDs matching Gene_ID in mapping CSV
- Genome FASTA headers must match GFF sequence names
"""
import os
import argparse
import logging
from collections import Counter
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

def get_mono_freq(seq):
    """Calculate mononucleotide frequencies f_X for A,C,G,T."""
    seq = seq.upper().replace("N", "")
    total = len(seq)
    if total == 0:
        return {b: 0.0 for b in "ACGT"}
    counts = Counter(seq)
    return {b: counts.get(b, 0) / total for b in "ACGT"}

def get_dinuc_counts(seq):
    """Count overlapping dinucleotides in sequence."""
    seq = seq.upper().replace("N", "")
    return Counter(seq[i:i+2] for i in range(len(seq) - 1) if len(seq[i:i+2]) == 2)

def karlin_rho(seq):
    """
    Compute Karlin's relative abundance vector rho for all 16 dinucleotides.
    rho_XY = f_XY / (f_X * f_Y)
    """
    mono = get_mono_freq(seq)
    dinuc_counts = get_dinuc_counts(seq)
    total_dinucs = sum(dinuc_counts.values())

    rho = {}
    for x in "ACGT":
        for y in "ACGT":
            dinuc = x + y
            f_xy = (dinuc_counts.get(dinuc, 0) / total_dinucs) if total_dinucs > 0 else 0.0
            f_x = mono.get(x, 0.0)
            f_y = mono.get(y, 0.0)
            if f_x > 0 and f_y > 0:
                rho_val = f_xy / (f_x * f_y)
            else:
                rho_val = 0.0
            rho[dinuc] = rho_val
    return rho

def karlin_delta_star(rho_gene, rho_genome):
    """Compute Karlin δ* as average absolute difference over 16 dinucleotides."""
    diffs = [abs(rho_gene[d] - rho_genome[d]) for d in sorted(rho_gene.keys())]
    return round(sum(diffs) / 16.0, 5)

def load_longest_contig_sequence(fna_file):
    """Load the longest contig sequence from genome FASTA."""
    longest_seq = ""
    max_len = 0
    for record in SeqIO.parse(fna_file, "fasta"):
        seq = str(record.seq)
        if len(seq) > max_len:
            longest_seq = seq
            max_len = len(seq)
    return longest_seq

def extract_gene_sequence_from_gff(genome_fna, gff_file, gene_id):
    """
    Extract gene CDS sequence by ID from GFF and genome FASTA.
    Returns sequence string, reverse complement if on '-' strand.
    """
    records = list(SeqIO.parse(genome_fna, "fasta"))
    chrom_dict = {rec.id: rec.seq for rec in records}

    with open(gff_file) as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            feature_type = parts[2].lower()
            if feature_type != "cds":
                continue
            attr_field = parts[8]
            if f"ID={gene_id}" not in attr_field:
                continue
            chrom = parts[0]
            start = int(parts[3]) - 1  # GFF is 1-based inclusive
            end = int(parts[4])        # inclusive
            strand = parts[6]
            contig_seq = chrom_dict.get(chrom)
            if contig_seq is None:
                logging.warning(f"Contig {chrom} not found in genome FASTA")
                return None
            subseq = contig_seq[start:end]
            if strand == "-":
                subseq = subseq.reverse_complement()
            return str(subseq).upper()
    return None

def compute_karlin_delta_for_genes(mapping_csv, genome_dir, output_csv, min_gene_length=30):
    """
    For each gene in mapping CSV, compute Karlin δ* against longest contig background.
    """
    df = pd.read_csv(mapping_csv)
    results = []
    genome_rho_cache = {}

    for idx, row in df.iterrows():
        genome = str(row['Genome'])
        gene_id = str(row['Gene_ID'])
        gene_name = str(row.get('Gene_Name', ""))

        fna_file = os.path.join(genome_dir, f"{genome}.fna")
        gff_file = os.path.join(genome_dir, f"{genome}.gff")

        if not os.path.isfile(fna_file) or not os.path.isfile(gff_file):
            logging.warning(f"Missing genome files for {genome}, skipping...")
            continue

        try:
            if genome not in genome_rho_cache:
                longest_contig_seq = load_longest_contig_sequence(fna_file)
                if len(longest_contig_seq) == 0:
                    logging.warning(f"No contigs found for genome {genome}")
                    continue
                genome_rho_cache[genome] = karlin_rho(longest_contig_seq)
            genome_rho = genome_rho_cache[genome]

            gene_seq = extract_gene_sequence_from_gff(fna_file, gff_file, gene_id)
            if gene_seq is None:
                logging.warning(f"Gene {gene_id} not found in genome {genome}")
                continue

            if len(gene_seq) < min_gene_length:
                logging.info(f"Skipping gene {gene_id} (length {len(gene_seq)}) — too short")
                continue

            gene_rho = karlin_rho(gene_seq)
            delta_star = karlin_delta_star(gene_rho, genome_rho)

            results.append({
                "Genome": genome,
                "Gene_ID": gene_id,
                "Gene_Name": gene_name,
                "Gene_Length": len(gene_seq),
                "Delta_Star": delta_star
            })

        except Exception as e:
            logging.error(f"Error processing gene {gene_id} in genome {genome}: {e}")

    out_df = pd.DataFrame(results)
    out_df.to_csv(output_csv, index=False)
    logging.info(f"Results saved to {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute Karlin δ* of genes vs longest contig background")
    parser.add_argument("-m", "--mapping", required=True, help="CSV file with Genome, Gene_ID, Gene_Name")
    parser.add_argument("-g", "--genome_dir", required=True, help="Directory with genome .fna and .gff files")
    parser.add_argument("-o", "--output", default="karlin_delta_star_results.csv", help="Output CSV file")
    parser.add_argument("-l", "--min_gene_length", type=int, default=30, help="Minimum gene length to analyze")
    args = parser.parse_args()

    compute_karlin_delta_for_genes(
        mapping_csv=args.mapping,
        genome_dir=args.genome_dir,
        output_csv=args.output,
        min_gene_length=args.min_gene_length
    )
