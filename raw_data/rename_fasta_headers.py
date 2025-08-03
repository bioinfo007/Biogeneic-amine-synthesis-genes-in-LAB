#!/usr/bin/env python3
"""
rename_headers.py

Usage:
    Make sure you have the following:
    - A CSV metadata file mapping accession to BioSample IDs (default: metadata_file.csv)
      CSV Format: accession,biosample_id (as column names)
    - A folder containing .fna files (default: HQ_genome)
    - An output folder to save renamed files (default: modified_HQ_genome)

Run the script:
    ./rename_headers.py
    OR
    python3 rename_headers.py

What it does:
    For each .fna file in the input folder, the script finds a matching accession
    in the CSV file. It then rewrites the FASTA headers by replacing them with
    "<BioSampleID>_<contig_number>" and writes the updated sequences to the output folder.

Notes:
    - Make the script executable with: chmod +x rename_headers.py
    - Modify input/output file/folder names inside the script as needed.
"""

import os
import csv
from pathlib import Path

def main():
    # Input parameters
    csv_file = "metadata_file.csv"
    fna_folder = "HQ_genome"  # Folder containing .fna files
    output_folder = "modified_HQ_genome"

    os.makedirs(output_folder, exist_ok=True)

    accession_to_biosample = {}
    with open(csv_file, newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) >= 2:
                accession, biosample = row[0].strip(), row[1].strip()
                accession_to_biosample[accession] = biosample

    for fna_file in Path(fna_folder).glob("*.fna"):
        fna_name = fna_file.stem
        matched_key = None

        for accession in accession_to_biosample:
            if accession in fna_name:
                matched_key = accession
                break

        if not matched_key:
            print(f"No match found for {fna_file.name}")
            continue

        biosample_id = accession_to_biosample[matched_key]
        contig_count = 0
        output_lines = []

        with open(fna_file, 'r') as infile:
            for line in infile:
                if line.startswith('>'):
                    contig_count += 1
                    new_header = f">{biosample_id}_{contig_count}"
                    output_lines.append(new_header + "\n")
                else:
                    output_lines.append(line)

        output_path = Path(output_folder) / fna_file.name
        with open(output_path, 'w') as outfile:
            outfile.writelines(output_lines)

        print(f"Renamed headers in {fna_file.name} → {output_path.name} ({contig_count} contigs)")

if __name__ == "__main__":
    main()
