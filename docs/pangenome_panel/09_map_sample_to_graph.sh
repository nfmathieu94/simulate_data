#!/usr/bin/bash -l
#SBATCH -p epyc
#SBATCH --job-name=pangenome_giraffe
#SBATCH --mem=64gb
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH -o logs/pangenome_panel/giraffe.%A_%a.log

# Map the Phase 2 sample reads onto the pangenome graph with vg giraffe.
# Array tasks: 0 Illumina (paired), 1 ONT-HQ, 2 PacBio HiFi.
#
# Giraffe handles both short and long reads as of 2025, so one mapper covers
# every technology in the panel.

set -euo pipefail

BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"

THREADS="${SLURM_CPUS_PER_TASK:-1}"
CONFIG="${CONFIG:-config/pangenome_panel.toml}"
OUTPUT="${OUTPUT:-results/pangenome_panel}"
GRAPH_DIR="$OUTPUT/graph"
GRAPH_NAME="${GRAPH_NAME:-riceTEpan}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

module load cactus/3.2.0    # provides vg 1.74.0

SAMPLE_NAME="$(
    python3.12 -c '
import sys, tomllib
with open(sys.argv[1], "rb") as handle:
    print(tomllib.load(handle).get("sample", {}).get("name", "SampleA"))
' "$CONFIG"
)"

SAMPLE_DIR="$OUTPUT/sample/$SAMPLE_NAME"
READS_ROOT="$SAMPLE_DIR/reads"

[[ -e "$GRAPH_DIR/.graph.complete" ]] || {
    echo "ERROR: pangenome graph is incomplete: $GRAPH_DIR" >&2
    exit 1
}

case "$TASK_ID" in
    0) TECH="illumina" ;;
    1) TECH="ont-hq" ;;
    2) TECH="hifi" ;;
    *) echo "ERROR: unknown array task: $TASK_ID" >&2; exit 1 ;;
esac

MAP_DIR="$SAMPLE_DIR/graph_alignments/$TECH"
SENTINEL="$MAP_DIR/.complete"
GAM="$MAP_DIR/${SAMPLE_NAME}.${TECH}.gam"

echo "[$(date)] Mapping $TECH reads to the pangenome graph"
echo "THREADS=$THREADS GRAPH=$GRAPH_NAME"

if [[ -e "$SENTINEL" ]]; then
    echo "[$(date)] $TECH already mapped; skipping"
    exit 0
fi
if [[ -d "$MAP_DIR" ]] && [[ -n "$(ls -A "$MAP_DIR" 2>/dev/null)" ]]; then
    echo "ERROR: refusing incomplete alignment directory: $MAP_DIR" >&2
    exit 1
fi
mkdir -p "$MAP_DIR"

GBZ="$GRAPH_DIR/${GRAPH_NAME}.gbz"
DIST="$GRAPH_DIR/${GRAPH_NAME}.dist"
# cactus-pangenome --giraffe emits a short-read minimizer index plus zipcodes.
MIN="$(ls "$GRAPH_DIR"/${GRAPH_NAME}*.min 2>/dev/null | head -1 || true)"
ZIP="$(ls "$GRAPH_DIR"/${GRAPH_NAME}*.zipcodes 2>/dev/null | head -1 || true)"

for required in "$GBZ" "$DIST"; do
    [[ -s "$required" ]] || {
        echo "ERROR: graph index is missing: $required" >&2
        exit 1
    }
done
[[ -n "$MIN" ]] || {
    echo "ERROR: no minimizer index found in $GRAPH_DIR" >&2
    exit 1
}

echo "GBZ=$GBZ"
echo "DIST=$DIST"
echo "MIN=$MIN"
echo "ZIP=${ZIP:-<none>}"

GIRAFFE_ARGS=(
    --gbz-name "$GBZ"
    --dist-name "$DIST"
    --minimizer-name "$MIN"
    --sample "$SAMPLE_NAME"
    --threads "$THREADS"
    --progress
)
[[ -n "$ZIP" ]] && GIRAFFE_ARGS+=(--zipcode-name "$ZIP")

case "$TECH" in
    illumina)
        R1="$READS_ROOT/illumina/illumina_reads1.fq.gz"
        R2="$READS_ROOT/illumina/illumina_reads2.fq.gz"
        for required in "$R1" "$R2"; do
            [[ -s "$required" ]] || {
                echo "ERROR: reads missing: $required" >&2
                exit 1
            }
        done
        vg giraffe "${GIRAFFE_ARGS[@]}" -f "$R1" -f "$R2" > "$GAM"
        ;;
    ont-hq|hifi)
        READS="$(ls "$READS_ROOT/$TECH"/*.f*q.gz 2>/dev/null | head -1 || true)"
        [[ -n "$READS" && -s "$READS" ]] || {
            echo "ERROR: no reads found under $READS_ROOT/$TECH" >&2
            exit 1
        }
        echo "READS=$READS"
        # Long reads need the long-read preset; short-read defaults produce
        # poor alignments for multi-kb reads.
        vg giraffe "${GIRAFFE_ARGS[@]}" --parameter-preset lr -f "$READS" > "$GAM"
        ;;
esac

[[ -s "$GAM" ]] || {
    echo "ERROR: giraffe produced an empty GAM: $GAM" >&2
    exit 1
}

echo
echo "[$(date)] Mapping complete; computing statistics"
vg stats -a "$GAM" > "$MAP_DIR/${SAMPLE_NAME}.${TECH}.stats.txt" 2>&1 || true
cat "$MAP_DIR/${SAMPLE_NAME}.${TECH}.stats.txt"

gzip -6 "$GAM"

touch "$SENTINEL"

echo
echo "[$(date)] $TECH graph mapping complete"
du -sh "$MAP_DIR"
