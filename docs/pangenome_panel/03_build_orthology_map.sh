#!/usr/bin/bash -l
#SBATCH -p epyc
#SBATCH --job-name=pangenome_orthology
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=16
#SBATCH --time=08:00:00
#SBATCH -o logs/pangenome_panel/orthology.%j.log

# Align the reference genome's panel chromosome against every other genome so
# anchor sites can be lifted between them.
#
# ARGUMENT ORDER IS LOAD-BEARING. minimap2 takes <target> <query>, and
# `paftools.js liftover` maps QUERY coordinates into TARGET coordinates.
# Anchors are chosen in the reference (Nipponbare), so the reference must be
# the QUERY and each other genome the TARGET. Reversing this silently
# produces a coordinate map that runs the wrong way.

set -euo pipefail

BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"

THREADS="${SLURM_CPUS_PER_TASK:-1}"
CONFIG="${CONFIG:-config/pangenome_panel.toml}"
OUTPUT="${OUTPUT:-results/pangenome_panel}"
CHR_DIR="$OUTPUT/genomes/chr1"
ORTHO_DIR="$OUTPUT/orthology"
STATS="$ORTHO_DIR/alignment_stats.tsv"

echo "[$(date)] Building orthology map"
echo "THREADS=$THREADS"

module load minimap2/2.30

mkdir -p "$ORTHO_DIR"

read -r REFERENCE PRESET MIN_LEN < <(
    python3.12 -c '
import sys, tomllib
with open(sys.argv[1], "rb") as handle:
    config = tomllib.load(handle)
reference = next(g["name"] for g in config["genomes"] if g.get("reference"))
ortho = config["orthology"]
print(reference, ortho["minimap2_preset"], ortho["min_alignment_length"])
' "$CONFIG"
)

mapfile -t TARGETS < <(
    python3.12 -c '
import sys, tomllib
with open(sys.argv[1], "rb") as handle:
    config = tomllib.load(handle)
for genome in config["genomes"]:
    if not genome.get("reference"):
        print(genome["name"])
' "$CONFIG"
)

QUERY="$CHR_DIR/${REFERENCE}.chr1.fa"
[[ -s "$QUERY" ]] || {
    echo "ERROR: reference chromosome is missing: $QUERY" >&2
    exit 1
}

REFERENCE_LENGTH="$(cut -f2 "${QUERY}.fai")"
printf 'query\ttarget\tblocks\taligned_bp\tlarge_blocks\tsyntenic_fraction\n' > "$STATS.tmp"

for TARGET in "${TARGETS[@]}"; do
    TARGET_FA="$CHR_DIR/${TARGET}.chr1.fa"
    PAF="$ORTHO_DIR/${REFERENCE}_to_${TARGET}.paf"
    SENTINEL="$ORTHO_DIR/.${TARGET}.complete"

    [[ -s "$TARGET_FA" ]] || {
        echo "ERROR: target chromosome is missing: $TARGET_FA" >&2
        exit 1
    }

    if [[ -e "$SENTINEL" && -s "$PAF" ]]; then
        echo "[$(date)] $REFERENCE -> $TARGET already complete; skipping"
    else
        echo "[$(date)] Aligning $REFERENCE (query) to $TARGET (target)"
        # -c is required, not optional: it performs base-level alignment and
        # emits the cg (CIGAR) tag. `paftools.js liftover` fails outright with
        # "unable to find the 'cg' tag" without it.
        minimap2 -c -x "$PRESET" -t "$THREADS" --cs \
            "$TARGET_FA" "$QUERY" > "$PAF.tmp"
        mv "$PAF.tmp" "$PAF"
        touch "$SENTINEL"
    fi

    # Column 11 is the alignment block length in PAF.
    awk -v q="$REFERENCE" -v t="$TARGET" -v reflen="$REFERENCE_LENGTH" \
        -v minlen="$MIN_LEN" '
        { blocks++; aligned += $11; if ($11 >= minlen) { large++; largebp += $11 } }
        END {
            printf "%s\t%s\t%d\t%d\t%d\t%.4f\n",
                q, t, blocks, aligned, large, (reflen ? largebp / reflen : 0)
        }
    ' "$PAF" >> "$STATS.tmp"
done

mv "$STATS.tmp" "$STATS"

echo
echo "[$(date)] Orthology map complete"
column -t -s $'\t' "$STATS"
echo
echo "Expect the highest syntenic fraction for the closest relative and the"
echo "lowest for the most diverged genome; the lowest one gates anchor yield."
