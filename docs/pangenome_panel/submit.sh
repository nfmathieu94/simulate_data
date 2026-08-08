#!/usr/bin/bash -l
# Submit Phase 1 of the multi-genome TE panel as a dependency chain:
#   download -> prepare chromosomes -> orthology -> annotation -> panel
# Run from the project root.

set -euo pipefail

CONFIG="${CONFIG:-config/pangenome_panel.toml}"
OUTPUT="${OUTPUT:-results/pangenome_panel}"
PIPELINE="pipeline/make_pangenome_panel"

[[ -f "$CONFIG" ]] || {
    echo "ERROR: config is missing: $CONFIG" >&2
    exit 1
}

mkdir -p logs/pangenome_panel

GENOME_COUNT="$(
    python3.12 -c '
import sys, tomllib
with open(sys.argv[1], "rb") as handle:
    print(len(tomllib.load(handle)["genomes"]))
' "$CONFIG"
)"
LAST_GENOME=$((GENOME_COUNT - 1))

export CONFIG OUTPUT

echo "Submitting download"
DOWNLOAD="$(sbatch --parsable "$PIPELINE/01_download_genomes.sh")"
echo "  $DOWNLOAD"

echo "Submitting chromosome preparation"
PREPARE="$(sbatch --parsable --dependency="afterok:${DOWNLOAD}" \
    "$PIPELINE/02_prepare_chromosomes.sh")"
echo "  $PREPARE"

echo "Submitting orthology map"
ORTHOLOGY="$(sbatch --parsable --dependency="afterok:${PREPARE}" \
    "$PIPELINE/03_build_orthology_map.sh")"
echo "  $ORTHOLOGY"

echo "Submitting background TE annotation (array 0-${LAST_GENOME})"
ANNOTATE="$(sbatch --parsable --dependency="afterok:${PREPARE}" \
    --array="0-${LAST_GENOME}" "$PIPELINE/04_annotate_background_tes.sh")"
echo "  $ANNOTATE"

echo "Submitting panel build"
PANEL="$(sbatch --parsable \
    --dependency="afterok:${ORTHOLOGY}:${ANNOTATE}" \
    "$PIPELINE/05_build_panel.sh")"
echo "  $PANEL"

echo
echo "Phase 1 submitted. Monitor with: squeue -u \"\$USER\""
echo "Phase 1 ends with the augmented genomes and truth table; no graph is built."
