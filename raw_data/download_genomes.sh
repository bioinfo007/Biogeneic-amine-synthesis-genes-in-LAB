#!/bin/bash
# Download all LAB genome assemblies using accession list, but keep only .fna and seq-report

set -euo pipefail

ACCESSION_FILE="accession_list.txt"
OUT_DIR="ncbi_genomes"

mkdir -p "$OUT_DIR"

while read -r accession; do
    echo "Processing $accession..."
    ZIP_PATH="${OUT_DIR}/${accession}.zip"
    TARGET_DIR="${OUT_DIR}/${accession}"

    # Download
    datasets download genome accession "$accession" --include genome,seq-report --filename "$ZIP_PATH"

    # Unzip (overwrite if exists)
    unzip -o "$ZIP_PATH" -d "$TARGET_DIR"
    rm "$ZIP_PATH"

    # Find and preserve .fna; move them to a clean subfolder, remove everything else
    PRESERVE_DIR="${TARGET_DIR}/kept"
    mkdir -p "$PRESERVE_DIR"

    # Move .fna files (could be multiple) and seq-report
    find "$TARGET_DIR" -type f \( -iname "*.fna" -o -iname "*seq-report*" \) -exec mv -t "$PRESERVE_DIR" {} +

    # Clean up everything except the preserved files
    find "$TARGET_DIR" -mindepth 1 -maxdepth 1 ! -name kept -exec rm -rf {} +

    # Flatten: move preserved files up one level and remove the extra folder
    mv "${PRESERVE_DIR}"/* "$TARGET_DIR"/
    rmdir "${PRESERVE_DIR}"

    echo "Retained files for $accession:"
    ls -1 "$TARGET_DIR"
done < "$ACCESSION_FILE"

    #move all genomes to HQ_genome folder
    mkdir -p HQ_genome
    mv ncbi_genomes/GC*/*.fna HQ_genome/
    rm -rf ncbi_genomes

echo "All downloads and pruning completed."
echo "Remove intermediate files and folders"
