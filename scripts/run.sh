#!/usr/bin/bash -l
set -euo pipefail

#SBATCH --job-name=sim_fastq
#SBATCH -p short
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=00:10:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1

BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"

mkdir -p logs results

pixi run simulate-data fastq \
    -n 1000 \
    -l 150 \
    --seed 42 \
    -o results/simulated_reads.fastq
