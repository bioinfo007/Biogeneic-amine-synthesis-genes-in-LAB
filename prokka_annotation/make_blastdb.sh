#!/usr/bin/env bash
set -euo pipefail

: <<'USAGE'
combine_and_make_blastdb.sh
---------------------------

Description:
    Given a path to a Prokka results root, recursively finds all .faa files under it,
    concatenates them into one FASTA (with sample-prefixed headers to avoid ID collisions), and builds
    a BLASTP database.

Usage:
    ./combine_and_make_blastdb.sh /path/to/prokka_results [combined_faa] [db_name]

    Arguments:
      1. Path to Prokka root directory (absolute or relative) (required)
      2. Combined .faa output filename (optional, default: combined_proteins.faa)
      3. BLASTP database name (optional, default: prokka_proteins_db)

Requirements:
    - BLAST+ (makeblastdb) in PATH

Outputs:
    - Combined FASTA file of all proteins (combined_faa)
    - BLASTP database files inside blastp_database/ directory with prefix db_name
USAGE

# Validate args
if [[ $# -lt 1 ]]; then
    echo "ERROR: Need path to Prokka results root." >&2
    echo "Usage: $0 /path/to/prokka_results [combined_faa] [db_name]" >&2
    exit 1
fi

PROKKA_ROOT_INPUT="$1"
COMBINED_FAA=${2:-combined_proteins.faa}
DB_NAME=${3:-prokka_proteins_db}

# Resolve path (works for relative or absolute)
if ! PROKKA_ROOT=$(realpath "$PROKKA_ROOT_INPUT" 2>/dev/null); then
    echo "ERROR: Failed to resolve path: $PROKKA_ROOT_INPUT" >&2
    exit 1
fi

# Check dependencies
if ! command -v makeblastdb &>/dev/null; then
    echo "ERROR: makeblastdb not found in PATH. Install BLAST+." >&2
    exit 1
fi

# Verify directory
if [[ ! -d "$PROKKA_ROOT" ]]; then
    echo "ERROR: Directory does not exist: $PROKKA_ROOT" >&2
    exit 1
fi

# Collect .faa files recursively
mapfile -t faa_files < <(find "$PROKKA_ROOT" -type f -name "*.faa")
if [[ ${#faa_files[@]} -eq 0 ]]; then
    echo "ERROR: No .faa files found under $PROKKA_ROOT" >&2
    exit 1
fi

echo "Found ${#faa_files[@]} .faa files under $PROKKA_ROOT"
echo "Combining into $COMBINED_FAA"

# Combine with sample-prefixed headers
: > "$COMBINED_FAA"  # truncate/create

for faa in "${faa_files[@]}"; do
    sample=$(basename "$(dirname "$faa")")
    while IFS= read -r line; do
        if [[ $line == '>'* ]]; then
            echo ">${sample}__${line#>}" >> "$COMBINED_FAA"
        else
            echo "$line" >> "$COMBINED_FAA"
        fi
    done < "$faa"
done

# Prepare database output folder
DB_FOLDER="blastp_database"
mkdir -p "$DB_FOLDER"
DB_PATH="$DB_FOLDER/$DB_NAME"

echo "Building BLASTP database ('$DB_PATH') from $COMBINED_FAA"
makeblastdb -in "$COMBINED_FAA" -dbtype prot -out "$DB_PATH" -parse_seqids

echo "Done. Database files stored in $DB_FOLDER with prefix: $DB_NAME"
