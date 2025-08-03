#!/usr/bin/env bash
set -euo pipefail

: <<'USAGE'
blastp_batch.sh
---------------

Description:
    Runs BLASTP for all .fasta query files in a folder against a specified BLASTP database.
    Outputs tabular format (outfmt 6) results as separate .tsv files inside BA_gene_hits directory.

Usage:
    ./blastp_batch.sh [query_folder] [blastdb_prefix]

Arguments:
    1. Folder containing query .fasta files (default: ref/)
    2. BLASTP database prefix (default: blastp_database/combined_LAB_DB)

Example:
    ./blastp_batch.sh
    ./blastp_batch.sh my_queries my_blast_db

Requirements:
    - blastp in PATH
USAGE

# Defaults
DEFAULT_QUERY_DIR="BA_gene_queries"
DEFAULT_BLAST_DB="blastp_database/combined_LAB_DB"

QUERY_DIR="${1:-$DEFAULT_QUERY_DIR}"
BLAST_DB="${2:-$DEFAULT_BLAST_DB}"
OUTPUT_DIR="BA_gene_hits"

if [[ ! -d "$QUERY_DIR" ]]; then
    echo "ERROR: Query folder not found: $QUERY_DIR" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

if ! command -v blastp &>/dev/null; then
    echo "ERROR: blastp not found in PATH." >&2
    exit 1
fi

shopt -s nullglob
query_files=("$QUERY_DIR"/*.fasta)
if [[ ${#query_files[@]} -eq 0 ]]; then
    echo "ERROR: No .fasta files found in $QUERY_DIR" >&2
    exit 1
fi

for f in "${query_files[@]}"; do
    base=$(basename "$f")
    out="$OUTPUT_DIR/${base}.tsv"
    echo "Running blastp: query=$f db=$BLAST_DB → output=$out"
    blastp -query "$f" -db "$BLAST_DB" -evalue 1e-5 -outfmt 6 -out "$out" -max_target_seqs 100000
done

echo "All queries processed. Results in $OUTPUT_DIR/"
