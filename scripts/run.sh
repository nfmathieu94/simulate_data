#!/usr/bin/env bash
set -euo pipefail

#SBATCH --job-name=sim_fastq
#SBATCH --output=../logs/%x_%j.out
#SBATCH --error=../logs/%x_%j.err
#SBATCH --time=00:10:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

simulate-data fastq \
    -n 1000 \
    -l 150 \
    --seed 42 \
    -o results/simulated_reads.fastq
