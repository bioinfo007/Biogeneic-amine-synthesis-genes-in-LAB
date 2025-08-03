#!/usr/bin/env bash
set -euo pipefail

: <<'USAGE'
bash_run_prokka.sh
------------------

Description:
    Runs Prokka on all .fna genome files in ../raw_data/modified_HQ_genome.
    Extracts locus tag from the first FASTA header and uses it (sanitized) as --locustag.
    Outputs go into prokka_results/<sample_name>.

Usage:
    ./bash_run_prokka.sh [--force]

    --force    Overwrite existing Prokka output folders. Without it, existing outputs are skipped.

Requirements:
    - prokka in PATH
    - Input .fna files in ../raw_data/modified_HQ_genome
    - BioPerl (Prokka dependency)

Outputs:
    - prokka_results/<sample_name>/ with annotated files
USAGE

# Check Prokka exists
if ! command -v prokka &>/dev/null; then
    echo "ERROR: prokka not found in PATH. Activate environment or install it." >&2
    exit 1
fi

# Configurable
INPUT_DIR="../raw_data/modified_HQ_genome"
OUT_ROOT="prokka_results"
CPUS=28
OVERWRITE=false

# Simple flag parse
if [[ "${1:-}" == "--force" ]]; then
    OVERWRITE=true
fi

mkdir -p "$OUT_ROOT"
shopt -s nullglob

ran=0
skipped=0
failed=0

for f in "$INPUT_DIR"/*.fna; do
    name=$(basename "$f" .fna)

    header_line=$(grep -m1 "^>" "$f" || true)
    if [[ -z "$header_line" ]]; then
        echo "WARNING: no FASTA header in $f; skipping." >&2
        ((skipped++))
        continue
    fi

    # Sanitize locus tag: alphanumeric + underscore, limit length
    sam_id=$(echo "${header_line#>}" | tr -cd '[:alnum:]_' | cut -c1-10)
    if [[ -z "$sam_id" ]]; then
        echo "WARNING: extracted empty locus tag from header in $f; using fallback '$name'." >&2
        sam_id="$name"
    fi

    outdir="$OUT_ROOT/$name"
    if [[ -d "$outdir" && -n "$(ls -A "$outdir")" && "$OVERWRITE" != "true" ]]; then
        echo "Skipping $name: output exists (use --force to overwrite)"
        ((skipped++))
        continue
    fi

    echo "Running Prokka on $f → $outdir (prefix=$name, locustag=$sam_id)"
    if prokka --outdir "$outdir" \
              --prefix "$name" \
              --locustag "$sam_id" \
              --cpus "$CPUS" \
              $( [[ "$OVERWRITE" == "true" ]] && echo "--force" ) \
              "$f"; then
        ((ran++))
    else
        echo "ERROR: Prokka failed for $name" >&2
        ((failed++))
    fi
done

# Summary
echo
echo "Summary:"
echo "  Annotated: $ran"
echo "  Skipped:   $skipped"
echo "  Failed:    $failed"
