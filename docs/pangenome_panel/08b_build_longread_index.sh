#!/usr/bin/bash -l
#SBATCH -p epyc
#SBATCH --job-name=pangenome_lr_index
#SBATCH --mem=64gb
#SBATCH --cpus-per-task=16
#SBATCH --time=08:00:00
#SBATCH -o logs/pangenome_panel/lr_index.%j.log

# Build a long-read minimizer index for the pangenome graph.
#
# cactus-pangenome --giraffe emits only a SHORT-READ minimizer index
# (*.shortread.withzip.min), built for ~150 bp reads. vg giraffe's long-read
# presets (hifi, r10) expect an index built with long-read (k, w) parameters,
# so multi-kb reads need their own index.
#
# vg also enforces that the minimizer index be NEWER than the distance index
# it depends on, and refuses to load an index that looks stale:
#   error[vg giraffe] <graph>.dist is newer than <graph>.min which depends on it
# Building here naturally satisfies that ordering.

set -euo pipefail

BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"

THREADS="${SLURM_CPUS_PER_TASK:-1}"
OUTPUT="${OUTPUT:-results/pangenome_panel}"
GRAPH_DIR="$OUTPUT/graph"
GRAPH_NAME="${GRAPH_NAME:-riceTEpan}"

# vg minimizer defaults (k=29, w=11) are the long-read settings; the
# short-read index cactus builds uses smaller kmers.
KMER="${KMER:-29}"
WINDOW="${WINDOW:-11}"

GBZ="$GRAPH_DIR/${GRAPH_NAME}.gbz"
DIST="$GRAPH_DIR/${GRAPH_NAME}.dist"
MIN="$GRAPH_DIR/${GRAPH_NAME}.longread.min"
ZIP="$GRAPH_DIR/${GRAPH_NAME}.longread.zipcodes"
SENTINEL="$GRAPH_DIR/.longread_index.complete"

echo "[$(date)] Building long-read minimizer index (k=$KMER w=$WINDOW)"

[[ -e "$GRAPH_DIR/.graph.complete" ]] || {
    echo "ERROR: pangenome graph is incomplete" >&2
    exit 1
}
for required in "$GBZ" "$DIST"; do
    [[ -s "$required" ]] || {
        echo "ERROR: missing graph index: $required" >&2
        exit 1
    }
done

if [[ -e "$SENTINEL" && -s "$MIN" ]]; then
    echo "[$(date)] Long-read index already built; skipping"
    exit 0
fi

module load cactus/3.2.0

vg minimizer \
    --distance-index "$DIST" \
    --output-name "$MIN" \
    --zipcode-name "$ZIP" \
    --kmer-length "$KMER" \
    --window-length "$WINDOW" \
    --threads "$THREADS" \
    --progress \
    "$GBZ"

[[ -s "$MIN" ]] || {
    echo "ERROR: vg minimizer produced no index: $MIN" >&2
    exit 1
}

# Guarantee the freshness ordering vg checks at load time.
touch "$MIN" "$ZIP"

touch "$SENTINEL"

echo
echo "[$(date)] Long-read index complete"
ls -la "$MIN" "$ZIP"
