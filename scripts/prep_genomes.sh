#!/usr/bin/env bash
# prep_genomes.sh — download benchmark reference genomes from NCBI and rename
# them to the labels used in benchmark/truth.tsv (so the harness can match them).
#
# Requires: NCBI `datasets` CLI  (conda install -c conda-forge ncbi-datasets-cli)
# Usage:    bash scripts/prep_genomes.sh benchmark/accessions.tsv benchmark/genomes
#
# Output:   <outdir>/<label>.fna   (nucleotide assembly per genome)
# NOTE:     verify accessions in accessions.tsv against NCBI before trusting —
#           some (esp. Leptospirillum) may need updating to a current assembly.
set -euo pipefail

ACC_TSV="${1:-benchmark/accessions.tsv}"
OUTDIR="${2:-benchmark/genomes}"
mkdir -p "$OUTDIR"

command -v datasets >/dev/null || { echo "[ERROR] NCBI 'datasets' CLI not found. conda install -c conda-forge ncbi-datasets-cli"; exit 1; }

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

grep -v '^#' "$ACC_TSV" | while IFS=$'\t' read -r label acc pheno; do
    [ -z "${label:-}" ] && continue
    out="$OUTDIR/${label}.fna"
    if [ -s "$out" ]; then echo "[skip] $label (exists)"; continue; fi
    echo "[get ] $label  $acc"
    datasets download genome accession "$acc" --include genome \
        --filename "$tmp/${label}.zip" >/dev/null
    unzip -oq "$tmp/${label}.zip" -d "$tmp/${label}"
    # concatenate all fasta parts into one labelled file
    find "$tmp/${label}" -name '*.fna' -exec cat {} + > "$out"
    [ -s "$out" ] || { echo "[WARN] empty download for $label ($acc)"; }
done

echo "[done] genomes in $OUTDIR/"
echo "Next: call ORFs for accuracy mode, e.g."
echo "  conda activate efesto"
echo "  python scripts/benchmark_fegenie.py --mode accuracy --fna_dir $OUTDIR \\"
echo "         --hmm_dir hmm_library --out benchmark/results --threads 8"
