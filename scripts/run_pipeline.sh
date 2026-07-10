#!/usr/bin/bash -l
set -euo pipefail

#SBATCH --job-name=sim_data
#SBATCH -p batch
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

# Example pipeline: TE insertion → SV placement → read simulation
# Adjust parameters and paths as needed.

BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"

mkdir -p logs results

REF="data/ref_genome/MSU_r7.fa"
TE_FA="data/TE_lib/mping.fa"
KNOWN_DEL="data/ref_genome/MSU_r7.fa.RepeatMasker.out"

for input in "$REF" "$TE_FA" "$KNOWN_DEL"; do
    if [[ ! -f "$input" ]]; then
        echo "Required input not found: $input" >&2
        exit 1
    fi
done

# Step 1: Insert 100 mPing TEs into Chr1-Chr5
pixi run simulate-data te-insertion \
    --ref "$REF" \
    --te "$TE_FA" \
    --known-del "$KNOWN_DEL" \
    --num 100 \
    --chroms "Chr1-Chr5" \
    --seed 42 \
    --output results/te_insertion/

# Step 2: Place 50 SVs in the TE-modified genome
pixi run simulate-data sv-placement \
    --ref results/te_insertion/final_genome.fa \
    --num-sv 50 \
    --sv-types DEL,DUP,INV,TRA \
    --seed 123 \
    --output results/sv_placement/

# Step 3: Simulate Illumina paired-end reads (20x coverage)
pixi run simulate-data reads-illumina \
    --ref results/sv_placement/final_genome.fa \
    --read-length 150 \
    --coverage 20 \
    --seed 456 \
    --output results/reads_illumina/

# Step 4: Simulate ONT long reads (20x coverage)
pixi run simulate-data reads-ont \
    --ref results/sv_placement/final_genome.fa \
    --coverage 20 \
    --seed 789 \
    --output results/reads_ont/

# Step 5: Simulate PacBio HiFi reads (20x coverage)
pixi run simulate-data reads-pacbio \
    --ref results/sv_placement/final_genome.fa \
    --coverage 20 \
    --read-type HiFi \
    --seed 101 \
    --output results/reads_pacbio/

echo "Pipeline complete. Results in results/"
