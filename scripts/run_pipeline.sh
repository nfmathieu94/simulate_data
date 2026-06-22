#!/usr/bin/env bash
set -euo pipefail

#SBATCH --job-name=sim_data
#SBATCH --output=../logs/%x_%j.out
#SBATCH --error=../logs/%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

# Example pipeline: TE insertion → SV placement → read simulation
# Adjust parameters and paths as needed.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

REF="data/ref_genome/MSU_r7.fa"
TE_FA="data/TE_lib/mping.fa"

# Step 1: Insert 100 mPing TEs into Chr1-Chr5
simulate-data te-insertion \
    --ref "$REF" \
    --te "$TE_FA" \
    --num 100 \
    --chroms "Chr1-Chr5" \
    --seed 42 \
    --output results/te_insertion/

# Step 2: Place 50 SVs in the TE-modified genome
simulate-data sv-placement \
    --ref results/te_insertion/modified_genome.fa \
    --num-sv 50 \
    --sv-types DEL,DUP,INV,TRA \
    --seed 123 \
    --output results/sv_placement/

# Step 3: Simulate Illumina paired-end reads (20x coverage)
simulate-data reads-illumina \
    --ref results/sv_placement/modified_genome.fa \
    --read-length 150 \
    --coverage 20 \
    --seed 456 \
    --output results/reads_illumina/

# Step 4: Simulate ONT long reads (20x coverage)
simulate-data reads-ont \
    --ref results/sv_placement/modified_genome.fa \
    --coverage 20 \
    --seed 789 \
    --output results/reads_ont/

# Step 5: Simulate PacBio HiFi reads (20x coverage)
simulate-data reads-pacbio \
    --ref results/sv_placement/modified_genome.fa \
    --coverage 20 \
    --read-type HiFi \
    --seed 101 \
    --output results/reads_pacbio/

echo "Pipeline complete. Results in results/"
