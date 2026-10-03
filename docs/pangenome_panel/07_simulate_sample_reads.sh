#!/usr/bin/bash -l
#SBATCH -p epyc
#SBATCH --job-name=pangenome_sample_reads
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH -o logs/pangenome_panel/sample_reads.%A_%a.log

# Simulate reads from the Phase 2 sample genome. Array tasks:
#   0  Illumina paired 150 bp
#   1  ONT high-accuracy
#   2  PacBio HiFi
#
# All FASTQ output is gzip compressed; no uncompressed FASTQ is retained.

set -euo pipefail

BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"

THREADS="${SLURM_CPUS_PER_TASK:-1}"
SIM_PROJECT_DIR="/rhome/nmath020/bigdata/github/github_tools/data_sim/simulate_data"
CONFIG="${CONFIG:-config/pangenome_panel.toml}"
OUTPUT="${OUTPUT:-results/pangenome_panel}"
COVERAGE="${COVERAGE:-30}"
SEED="${SEED:-1017}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"

SAMPLE_NAME="$(
    python3.12 -c '
import sys, tomllib
with open(sys.argv[1], "rb") as handle:
    print(tomllib.load(handle).get("sample", {}).get("name", "SampleA"))
' "$CONFIG"
)"

SAMPLE_DIR="$OUTPUT/sample/$SAMPLE_NAME"
GENOME="$SAMPLE_DIR/${SAMPLE_NAME}.chr1.fa"

[[ -e "$SAMPLE_DIR/.complete" ]] || {
    echo "ERROR: sample genome is incomplete: $SAMPLE_DIR" >&2
    exit 1
}
[[ -s "$GENOME" ]] || {
    echo "ERROR: sample genome FASTA is missing: $GENOME" >&2
    exit 1
}

case "$TASK_ID" in
    0) TECH="illumina" ;;
    1) TECH="ont-hq" ;;
    2) TECH="hifi" ;;
    *) echo "ERROR: unknown array task: $TASK_ID" >&2; exit 1 ;;
esac

READ_DIR="$SAMPLE_DIR/reads/$TECH"
SENTINEL="$READ_DIR/.complete"

echo "[$(date)] Simulating $TECH reads for $SAMPLE_NAME"
echo "COVERAGE=${COVERAGE}x THREADS=$THREADS"

if [[ -e "$SENTINEL" ]]; then
    echo "[$(date)] $TECH reads already complete; skipping"
    exit 0
fi
if [[ -d "$READ_DIR" ]] && [[ -n "$(ls -A "$READ_DIR" 2>/dev/null)" ]]; then
    echo "ERROR: refusing incomplete read directory: $READ_DIR" >&2
    exit 1
fi
mkdir -p "$READ_DIR"

case "$TECH" in
    illumina)
        pixi run --manifest-path "$SIM_PROJECT_DIR/pyproject.toml" \
            simulate-data reads-illumina \
            --ref "$GENOME" \
            --coverage "$COVERAGE" \
            --read-length 150 \
            --fragment-size 300 \
            --fragment-std 30 \
            --sequencing-system HSXn \
            --seed "$SEED" \
            --output "$READ_DIR"
        ;;
    ont-hq)
        pixi run --manifest-path "$SIM_PROJECT_DIR/pyproject.toml" \
            simulate-data reads-ont \
            --ref "$GENOME" \
            --coverage "$COVERAGE" \
            --read-length 12000 \
            --read-std 9000 \
            --error-model QSHMM-ONT-HQ \
            --seed "$SEED" \
            --output "$READ_DIR"
        ;;
    hifi)
        pixi run --manifest-path "$SIM_PROJECT_DIR/pyproject.toml" \
            simulate-data reads-pacbio \
            --ref "$GENOME" \
            --coverage "$COVERAGE" \
            --read-type HiFi \
            --read-length 15000 \
            --read-std 2000 \
            --pass-num 10 \
            --threads "$THREADS" \
            --min-ccs-yield 0.5 \
            --seed "$SEED" \
            --output "$READ_DIR"
        ;;
esac

# ART writes uncompressed FASTQ; PBSIM3 and ccs already write gzip. Compress
# anything left so the whole panel keeps one contract.
shopt -s nullglob
for fastq in "$READ_DIR"/*.fq "$READ_DIR"/*.fastq; do
    gzip -6 "$fastq"
done
# ART's per-read alignment files are large and not part of the contract.
rm -f "$READ_DIR"/*.aln "$READ_DIR"/*.sam

# PBSIM3 MAF records every subread pass, so for HiFi it dwarfs the reads
# themselves (8.8 GB vs 1.2 GB at 30x). Nothing downstream consumes it here:
# sample truth is known by construction from sample_truth.tsv. Set
# KEEP_ALIGNMENTS=1 to retain MAF/ref for read-level provenance work.
if [[ "${KEEP_ALIGNMENTS:-0}" != "1" ]]; then
    rm -f "$READ_DIR"/*.maf.gz "$READ_DIR"/*.ref
fi
shopt -u nullglob

if find "$READ_DIR" -type f \( -name '*.fastq' -o -name '*.fq' \) \
    -print -quit | grep -q .; then
    echo "ERROR: uncompressed FASTQs remain in $READ_DIR" >&2
    exit 1
fi
if ! find "$READ_DIR" -type f -name '*.gz' -print -quit | grep -q .; then
    echo "ERROR: no compressed FASTQ produced in $READ_DIR" >&2
    exit 1
fi

touch "$SENTINEL"

echo
echo "[$(date)] $TECH read simulation complete"
du -sh "$READ_DIR"
ls -la "$READ_DIR"
