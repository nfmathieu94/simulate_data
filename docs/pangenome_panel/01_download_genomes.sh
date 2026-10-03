#!/usr/bin/bash -l
#SBATCH -p short
#SBATCH --job-name=pangenome_download
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH -o logs/pangenome_panel/download.%j.log

# Download the four MAGIC16 rice assemblies that make up the multi-genome TE
# panel, verifying each against NCBI's published MD5. Safe to rerun: a genome
# whose sentinel exists and whose MD5 still matches is skipped.

set -euo pipefail

BASE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$BASE_DIR"

CONFIG="${CONFIG:-config/pangenome_panel.toml}"
OUTPUT="${OUTPUT:-results/pangenome_panel}"
RAW_DIR="$OUTPUT/genomes/raw"
MANIFEST="$OUTPUT/genomes/download_manifest.tsv"

echo "[$(date)] Downloading rice assemblies for the pangenome panel"
echo "CONFIG=$CONFIG"
echo "OUTPUT=$OUTPUT"

[[ -f "$CONFIG" ]] || {
    echo "ERROR: config is missing: $CONFIG" >&2
    exit 1
}
command -v curl >/dev/null 2>&1 || {
    echo "ERROR: curl is unavailable" >&2
    exit 127
}

mkdir -p "$RAW_DIR"

# name<TAB>accession<TAB>directory, one genome per line
mapfile -t GENOMES < <(
    python3.12 -c '
import sys, tomllib
with open(sys.argv[1], "rb") as handle:
    config = tomllib.load(handle)
for genome in config["genomes"]:
    print("\t".join((genome["name"], genome["accession"], genome["directory"])))
' "$CONFIG"
)

BASE_URL="$(
    python3.12 -c '
import sys, tomllib
with open(sys.argv[1], "rb") as handle:
    print(tomllib.load(handle)["download"]["base_url"])
' "$CONFIG"
)"

printf 'name\taccession\turl\tbytes\tmd5\tdownloaded_at\n' > "$MANIFEST.tmp"

for record in "${GENOMES[@]}"; do
    NAME="$(cut -f1 <<<"$record")"
    ACCESSION="$(cut -f2 <<<"$record")"
    DIRECTORY="$(cut -f3 <<<"$record")"

    # NCBI splits the numeric part of the accession into three-digit path
    # components: GCA_009830595.1 -> GCA/009/830/595
    DIGITS="${ACCESSION#GCA_}"
    DIGITS="${DIGITS%%.*}"
    URL_DIR="$BASE_URL/${DIGITS:0:3}/${DIGITS:3:3}/${DIGITS:6:3}/$DIRECTORY"

    GENOME_DIR="$RAW_DIR/$NAME"
    FASTA="$GENOME_DIR/${DIRECTORY}_genomic.fna.gz"
    REPORT="$GENOME_DIR/${DIRECTORY}_assembly_report.txt"
    CHECKSUMS="$GENOME_DIR/md5checksums.txt"
    SENTINEL="$GENOME_DIR/.download.complete"

    mkdir -p "$GENOME_DIR"

    echo "[$(date)] $NAME ($ACCESSION)"

    if [[ -e "$SENTINEL" && -s "$FASTA" ]]; then
        EXPECTED="$(awk -v f="./${DIRECTORY}_genomic.fna.gz" '$2 == f {print $1}' \
            "$CHECKSUMS" 2>/dev/null || true)"
        OBSERVED="$(md5sum "$FASTA" | cut -d' ' -f1)"
        if [[ -n "$EXPECTED" && "$EXPECTED" == "$OBSERVED" ]]; then
            echo "  already complete, MD5 verified; skipping"
            printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$NAME" "$ACCESSION" "$URL_DIR" \
                "$(stat -c %s "$FASTA")" "$OBSERVED" "cached" >> "$MANIFEST.tmp"
            continue
        fi
        echo "  sentinel present but MD5 does not match; re-downloading" >&2
        rm -f "$SENTINEL"
    fi

    # --continue-at resumes a partial transfer instead of truncating it.
    curl --fail --location --silent --show-error \
        --output "$CHECKSUMS" "$URL_DIR/md5checksums.txt"
    curl --fail --location --silent --show-error \
        --output "$REPORT" "$URL_DIR/${DIRECTORY}_assembly_report.txt"
    curl --fail --location --show-error --continue-at - \
        --output "$FASTA" "$URL_DIR/${DIRECTORY}_genomic.fna.gz"

    EXPECTED="$(awk -v f="./${DIRECTORY}_genomic.fna.gz" '$2 == f {print $1}' "$CHECKSUMS")"
    [[ -n "$EXPECTED" ]] || {
        echo "ERROR: no MD5 published for $NAME in $CHECKSUMS" >&2
        exit 1
    }
    OBSERVED="$(md5sum "$FASTA" | cut -d' ' -f1)"
    if [[ "$EXPECTED" != "$OBSERVED" ]]; then
        echo "ERROR: MD5 mismatch for $NAME" >&2
        echo "  expected $EXPECTED" >&2
        echo "  observed $OBSERVED" >&2
        exit 1
    fi
    echo "  MD5 verified: $OBSERVED"

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$NAME" "$ACCESSION" "$URL_DIR" \
        "$(stat -c %s "$FASTA")" "$OBSERVED" "$(date -Is)" >> "$MANIFEST.tmp"
    touch "$SENTINEL"
done

mv "$MANIFEST.tmp" "$MANIFEST"

echo
echo "[$(date)] Download complete"
column -t -s $'\t' "$MANIFEST"
