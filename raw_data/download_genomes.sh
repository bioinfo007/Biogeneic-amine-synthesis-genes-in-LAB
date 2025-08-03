#!/bin/bash
# Download all LAB genome assemblies from accession list

# Exit on error
set -euo pipefail

ACCESSION_FILE="accession_list.txt"
OUT_DIR="ncbi_genomes"

mkdir -p "$OUT_DIR"

# Loop through each accession and download
while read -r accession; do
    echo "Downloading $accession..."
    datasets download genome accession "$accession" --include genome,gff3,protein --filename "${OUT_DIR}/${accession}.zip"
    unzip -o "${OUT_DIR}/${accession}.zip" -d "${OUT_DIR}/${accession}"
    rm "${OUT_DIR}/${accession}.zip"
done < "$ACCESSION_FILE"

echo "All downloads completed."
