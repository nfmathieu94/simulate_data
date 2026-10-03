#!/usr/bin/bash -l
#SBATCH -p short
#SBATCH --job-name=pangenome_prepare
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH -o logs/pangenome_panel/prepare.%j.log

# Extract the panel chromosome from each assembly and rename it to PanSN form
# (<sample>#1#Chr1), which is what pangenome graph tools expect.

set -euo pipefail

BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"

SIM_PROJECT_DIR="/rhome/nmath020/bigdata/github/github_tools/data_sim/simulate_data"
CONFIG="${CONFIG:-config/pangenome_panel.toml}"
OUTPUT="${OUTPUT:-results/pangenome_panel}"
SCRIPT="pipeline/make_pangenome_panel/prepare_chromosomes.py"

echo "[$(date)] Extracting panel chromosomes and applying PanSN naming"
echo "CONFIG=$CONFIG"
echo "OUTPUT=$OUTPUT"

command -v pixi >/dev/null 2>&1 || {
    echo "ERROR: pixi is unavailable" >&2
    exit 127
}

pixi run --manifest-path "$SIM_PROJECT_DIR/pyproject.toml" python "$SCRIPT" \
    --config "$CONFIG" \
    --output "$OUTPUT"

module load samtools/1.22.1
CHR_DIR="$OUTPUT/genomes/chr1"
for fasta in "$CHR_DIR"/*.chr1.fa; do
    [[ -s "$fasta" ]] || {
        echo "ERROR: empty chromosome FASTA: $fasta" >&2
        exit 1
    }
    samtools faidx "$fasta"
done

echo
echo "[$(date)] Chromosome preparation complete"
cut -f1,2 "$CHR_DIR"/*.fai
