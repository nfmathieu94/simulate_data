#!/usr/bin/bash -l
#SBATCH -p epyc
#SBATCH --job-name=pangenome_graph
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH -o logs/pangenome_panel/graph.%j.log

# Build a pangenome graph from the four TE-augmented genomes with
# Minigraph-Cactus, and emit the Giraffe indexes needed to map sample reads.
#
# Minigraph-Cactus rather than PGGB: the 2025 GigaScience crop-plant
# evaluation could not even build mapping indexes for a PGGB graph of six
# sorghum assemblies due to memory, while Minigraph-Cactus + Giraffe held
# >98% correct mapping. minigraph alone would miss sub-50 bp variation, which
# includes most TSDs.

set -euo pipefail

BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"

THREADS="${SLURM_CPUS_PER_TASK:-1}"
CONFIG="${CONFIG:-config/pangenome_panel.toml}"
OUTPUT="${OUTPUT:-results/pangenome_panel}"
GRAPH_DIR="$OUTPUT/graph"
GRAPH_NAME="${GRAPH_NAME:-riceTEpan}"
SENTINEL="$GRAPH_DIR/.graph.complete"

# Cactus writes a large jobstore and scratch; keep both off /bigdata.
SCRATCH_ROOT="${SCRATCH:-/scratch/$USER}/cactus.${SLURM_JOB_ID:-$$}"
JOBSTORE="$SCRATCH_ROOT/jobstore"
WORK_DIR="$SCRATCH_ROOT/work"

echo "[$(date)] Building the pangenome graph"
echo "THREADS=$THREADS"
echo "GRAPH_DIR=$GRAPH_DIR"

[[ -e "$OUTPUT/.panel.complete" ]] || {
    echo "ERROR: panel is incomplete; run the Phase 1 build first" >&2
    exit 1
}

if [[ -e "$SENTINEL" ]]; then
    echo "[$(date)] Graph already complete; skipping"
    exit 0
fi
if [[ -d "$GRAPH_DIR" ]] && [[ -n "$(ls -A "$GRAPH_DIR" 2>/dev/null)" ]]; then
    echo "ERROR: refusing incomplete graph directory: $GRAPH_DIR" >&2
    echo "  remove it deliberately to rebuild" >&2
    exit 1
fi

module load cactus/3.2.0

mkdir -p "$GRAPH_DIR" "$WORK_DIR"

REFERENCE="$(
    python3.12 -c '
import sys, tomllib
with open(sys.argv[1], "rb") as handle:
    config = tomllib.load(handle)
print(next(g["name"] for g in config["genomes"] if g.get("reference")))
' "$CONFIG"
)"

# Cactus seqFile: <sampleName><TAB><path>, one genome per line. The FASTA
# records are already PanSN-named (<sample>#1#Chr1) with prefixes matching
# these sample names, which is what the graph tools expect.
SEQFILE="$GRAPH_DIR/${GRAPH_NAME}.seqfile"
python3.12 - "$CONFIG" "$OUTPUT" "$SEQFILE" <<'PYSEQ'
import os
import sys
import tomllib

config_path, output, seqfile = sys.argv[1:4]
with open(config_path, "rb") as handle:
    config = tomllib.load(handle)
chrom = config["panel"]["chrom"].lower()

with open(seqfile, "w") as out:
    for genome in config["genomes"]:
        name = genome["name"]
        path = os.path.abspath(
            os.path.join(output, "genomes", "augmented", f"{name}.{chrom}.te.fa")
        )
        if not os.path.isfile(path):
            raise SystemExit(f"ERROR: augmented genome missing: {path}")
        out.write(f"{name}\t{path}\n")
PYSEQ

echo "--- seqFile ---"
cat "$SEQFILE"
echo "REFERENCE=$REFERENCE"

# single_machine keeps the whole workflow inside this allocation rather than
# having Toil submit its own SLURM jobs.
cactus-pangenome "$JOBSTORE" "$SEQFILE" \
    --outDir "$GRAPH_DIR" \
    --outName "$GRAPH_NAME" \
    --reference "$REFERENCE" \
    --batchSystem single_machine \
    --maxCores "$THREADS" \
    --workDir "$WORK_DIR" \
    --gfa clip full \
    --gbz clip \
    --vcf clip \
    --giraffe clip \
    --logInfo

rm -rf "$SCRATCH_ROOT"

for required in \
    "$GRAPH_DIR/${GRAPH_NAME}.gbz" \
    "$GRAPH_DIR/${GRAPH_NAME}.dist"; do
    [[ -s "$required" ]] || {
        echo "ERROR: expected graph output is missing: $required" >&2
        exit 1
    }
done

touch "$SENTINEL"

echo
echo "[$(date)] Pangenome graph complete"
ls -la "$GRAPH_DIR"
echo
echo "--- graph statistics ---"
vg stats -z -l "$GRAPH_DIR/${GRAPH_NAME}.gbz" || true
