#!/usr/bin/bash -l
#SBATCH -p epyc
#SBATCH --job-name=pangenome_annotate
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH -o logs/pangenome_panel/annotate.%A_%a.log

# Annotate pre-existing riceTElib-family TEs in each genome's panel
# chromosome. Synthetic insertions are then kept clear of them, so a correct
# RelocaTE3 call is never scored against a real element of the same family.
#
# RepeatMasker against riceTElib is used rather than EDTA/panEDTA: we only
# need to know where riceTElib FAMILIES already sit, not to discover novel
# TEs, and de-novo annotation of four genomes would dominate Phase 1 runtime.

set -euo pipefail

BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"

THREADS="${SLURM_CPUS_PER_TASK:-1}"
CONFIG="${CONFIG:-config/pangenome_panel.toml}"
OUTPUT="${OUTPUT:-results/pangenome_panel}"
CHR_DIR="$OUTPUT/genomes/chr1"
ANNOT_DIR="$OUTPUT/annotation"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

module load RepeatMasker/4.1.8

mapfile -t GENOMES < <(
    python3.12 -c '
import sys, tomllib
with open(sys.argv[1], "rb") as handle:
    for genome in tomllib.load(handle)["genomes"]:
        print(genome["name"])
' "$CONFIG"
)
TE_LIBRARY="$(
    python3.12 -c '
import sys, tomllib
with open(sys.argv[1], "rb") as handle:
    print(tomllib.load(handle)["inputs"]["te_library"])
' "$CONFIG"
)"

GENOME="${GENOMES[$TASK_ID]}"
FASTA="$CHR_DIR/${GENOME}.chr1.fa"
GENOME_DIR="$ANNOT_DIR/$GENOME"
BED="$ANNOT_DIR/${GENOME}.background_te.bed"
SENTINEL="$ANNOT_DIR/.${GENOME}.complete"

echo "[$(date)] Annotating background TEs"
echo "GENOME=$GENOME THREADS=$THREADS"
echo "TE_LIBRARY=$TE_LIBRARY"

[[ -s "$FASTA" ]] || {
    echo "ERROR: chromosome FASTA is missing: $FASTA" >&2
    exit 1
}
[[ -s "$TE_LIBRARY" ]] || {
    echo "ERROR: TE library is missing: $TE_LIBRARY" >&2
    exit 1
}

if [[ -e "$SENTINEL" && -s "$BED" ]]; then
    echo "[$(date)] $GENOME already annotated; skipping"
    exit 0
fi

mkdir -p "$GENOME_DIR"

# RepeatMasker's -pa counts JOBS OF 4 THREADS, so divide the allocation.
PA=$(( THREADS / 4 ))
[[ "$PA" -ge 1 ]] || PA=1

# -engine rmblast is required: the module defaults to HMMER, which indexes
# Dfam HMM profiles and cannot hmmpress a nucleotide library like riceTElib.
# rmblastn is also the engine GraffiTE uses for the same job.
#
# RepeatMasker drops RM_* working directories into the current directory, so
# run it from scratch space rather than littering the project root.
WORK_DIR="${SCRATCH:-/scratch/$USER}/repeatmasker.${SLURM_JOB_ID:-$$}.${GENOME}"
mkdir -p "$WORK_DIR"
ABS_LIBRARY="$(readlink -f "$TE_LIBRARY")"
ABS_FASTA="$(readlink -f "$FASTA")"
ABS_OUTDIR="$(readlink -f "$GENOME_DIR")"

(
    cd "$WORK_DIR"
    RepeatMasker -engine rmblast -pa "$PA" -lib "$ABS_LIBRARY" \
        -dir "$ABS_OUTDIR" -nolow -no_is -gff "$ABS_FASTA"
)
rm -rf "$WORK_DIR"

OUT_FILE="$GENOME_DIR/$(basename "$FASTA").out"
[[ -s "$OUT_FILE" ]] || {
    echo "ERROR: RepeatMasker produced no .out file: $OUT_FILE" >&2
    exit 1
}

# RepeatMasker .out: 3 header lines, then whitespace-delimited records with
# query start/end in columns 6/7 and the matching repeat name in column 10.
awk 'NR > 3 && NF >= 10 {
    printf "%s\t%d\t%s\t%s\n", $5, $6 - 1, $7, $10
}' "$OUT_FILE" | sort -k1,1 -k2,2n > "$BED.tmp"
mv "$BED.tmp" "$BED"

touch "$SENTINEL"

echo "[$(date)] Background annotation complete: $GENOME"
echo "intervals: $(wc -l < "$BED")"
