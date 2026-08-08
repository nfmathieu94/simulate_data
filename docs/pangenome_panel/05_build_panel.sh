#!/usr/bin/bash -l
#SBATCH -p epyc
#SBATCH --job-name=pangenome_panel
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH -o logs/pangenome_panel/panel.%j.log

# Select orthologous anchors, assign them to sharing patterns, and write the
# TE-augmented genomes plus the truth table. This is the end of Phase 1.

set -euo pipefail

BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"

SIM_PROJECT_DIR="/rhome/nmath020/bigdata/github/github_tools/data_sim/simulate_data"
CONFIG="${CONFIG:-config/pangenome_panel.toml}"
OUTPUT="${OUTPUT:-results/pangenome_panel}"
BUILDER="pipeline/make_pangenome_panel/05_build_pangenome_panel.py"

echo "[$(date)] Building the multi-genome TE panel"
echo "CONFIG=$CONFIG OUTPUT=$OUTPUT"

module load minimap2/2.30    # provides paftools.js and k8 for liftover

command -v paftools.js >/dev/null 2>&1 || {
    echo "ERROR: paftools.js is unavailable" >&2
    exit 127
}

pixi run --manifest-path "$SIM_PROJECT_DIR/pyproject.toml" python "$BUILDER" \
    --config "$CONFIG" \
    --output "$OUTPUT"

[[ -e "$OUTPUT/.panel.complete" ]] || {
    echo "ERROR: panel completion sentinel is missing" >&2
    exit 1
}
for required in "$OUTPUT/truth_events.tsv" "$OUTPUT/sharing_matrix.tsv"; do
    [[ -s "$required" ]] || {
        echo "ERROR: missing output: $required" >&2
        exit 1
    }
done

echo
echo "[$(date)] Panel complete"
echo "events: $(( $(wc -l < "$OUTPUT/truth_events.tsv") - 1 ))"
