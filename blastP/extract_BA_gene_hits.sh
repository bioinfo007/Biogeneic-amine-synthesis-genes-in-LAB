#!/usr/bin/env bash
set -euo pipefail

: <<'USAGE'
extract_fasta_from_blast_hits.sh
--------------------------------

Description:
    Extract protein sequences from combined_LAB.faa based on locus tags
    found in BLASTP .tsv files under BA_gene_hits directory.

Usage:
    ./extract_fasta_from_blast_hits.sh [tsv_folder] [combined_faa] [output_folder]

Defaults (if no arguments):
    tsv_folder: BA_gene_hits
    combined_faa: combined_LAB.faa
    output_folder: extracted_BA_genes

Requirements:
    - Python 3 with Biopython installed
USAGE

TSV_DIR="${1:-BA_gene_hits}"
COMBINED_FASTA="${2:-combined_LAB.faa}"
OUT_DIR="${3:-extracted_BA_genes}"

if [[ ! -d "$TSV_DIR" ]]; then
    echo "ERROR: TSV directory not found: $TSV_DIR" >&2
    exit 1
fi

if [[ ! -f "$COMBINED_FASTA" ]]; then
    echo "ERROR: Combined FASTA file not found: $COMBINED_FASTA" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

extract_script=$(cat << 'EOF'
import sys
from Bio import SeqIO

tsv_file = sys.argv[1]
fasta_file = sys.argv[2]
output_file = sys.argv[3]

ids = set()
with open(tsv_file) as f:
    for line in f:
        if line.strip():
            parts = line.strip().split('\t')
            if len(parts) > 1:
                ids.add(parts[1])

seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))

with open(output_file, 'w') as out_f:
    for seqid in ids:
        if seqid in seq_dict:
            SeqIO.write(seq_dict[seqid], out_f, 'fasta')
EOF
)

shopt -s nullglob
tsv_files=("$TSV_DIR"/*.tsv)

if [[ ${#tsv_files[@]} -eq 0 ]]; then
    echo "ERROR: No .tsv files found in $TSV_DIR" >&2
    exit 1
fi

for tsv in "${tsv_files[@]}"; do
    base=$(basename "$tsv" .tsv)
    out_fasta="$OUT_DIR/${base}_hits.faa"
    echo "Extracting sequences for $base from $COMBINED_FASTA → $out_fasta"
    python3 -c "$extract_script" "$tsv" "$COMBINED_FASTA" "$out_fasta"
done

echo "Extraction done. Results saved in $OUT_DIR"
