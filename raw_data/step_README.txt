README
Overview
Two scripts to download genome files and rename FASTA headers based on metadata to be used for consistent naming in prokka.

Input Files
accession_list.txt: List of genome accession IDs (one per line).

metadata_file.csv: CSV mapping accession to BioSample ID, format: accession,biosample_id.

Scripts & Usage
1. download_genomes.sh
Downloads genomes using accession IDs from accession_list.txt.

Usage:
bash download_genomes.sh
Downloads saved to folder (e.g., HQ_genome/).

2. rename_fasta_headers.py
Renames headers in .fna files using BioSample IDs from metadata_file.csv.

Usage:
chmod +x rename_fasta_headers.py
./rename_fasta_headers.py
Outputs renamed files to modified_HQ_genome/.

Notes
Ensure accession IDs match across all files and filenames.
Ensure NCBI Datasets CLI is installed
Modify script variables if needed.
