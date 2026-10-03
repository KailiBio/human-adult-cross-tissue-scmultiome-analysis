#!/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE' >&2
Usage: generate_phylop_aggregation.sh -a CCRE_BED -w BIGWIG -o OUTFILE
                                       [-n N_BINS] [--window BP] [--column COL]
                                       [--chrom-sizes FILE] [--label LABEL]
                                       [--keep-matrix]

Generate a phyloP/conservation aggregation profile for a given list of
cCREs: recenters each cCRE to a fixed-width window (clipped so it doesn't
start before position 0), splits that window into N_BINS equal bins,
computes the mean bigWig signal (e.g. phyloP, phastCons, GERP) in each bin
per cCRE (bigWigAverageOverBed), and reduces across all cCREs to a single
mean aggregation profile ready for plotting.

Re-implementation (bash/awk/UCSC tools) of
0-7-2c_peak_analysis.conservation.sh (section 1's per-cCRE recentering,
daily/make_aggregation_matrix.sh's binning/signal/matrix steps, and
section 3's mean-profile reduction), generalized to any single cCRE list
and bigWig track. The project-specific overlap/tier stratifications and
random-region baseline are left to the caller: run this script again on a
filtered cCRE subset, or on a random-region BED, and compare the resulting
profiles downstream.

Requires bigWigAverageOverBed on PATH.

  -a CCRE_BED      cCRE BED, >=4 columns (chrom, start, end, id); id must
                   be unique
  -w BIGWIG        conservation (or other) bigWig track to aggregate
  -o OUTFILE       output path for the mean aggregation profile
  -n N_BINS        number of bins the window is split into (default: 400)
  --window BP      total window width, centered on each cCRE's midpoint
                   (default: 2000)
  --column COL     bigWigAverageOverBed output column to use as the
                   per-bin signal (default: 5, mean0 -- treats
                   non-covered bases as zero)
  --chrom-sizes    genome chrom.sizes; if given, windows are also clipped
                   so they don't run past a chromosome's end (the
                   original script only clipped the start)
  --label LABEL    row label for the output profile (default: basename of
                   -a without its extension)
  --keep-matrix    also keep the per-cCRE x bin signal matrix
                   (OUTFILE.matrix.txt), not just its mean
USAGE
  exit 1
}

N_BINS=400
WINDOW=2000
COLUMN=5
CHROM_SIZES=""
LABEL=""
KEEP_MATRIX=0
CCRE_BED=""
BIGWIG=""
OUTFILE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -a) CCRE_BED="$2"; shift 2 ;;
    -w) BIGWIG="$2"; shift 2 ;;
    -o) OUTFILE="$2"; shift 2 ;;
    -n) N_BINS="$2"; shift 2 ;;
    --window) WINDOW="$2"; shift 2 ;;
    --column) COLUMN="$2"; shift 2 ;;
    --chrom-sizes) CHROM_SIZES="$2"; shift 2 ;;
    --label) LABEL="$2"; shift 2 ;;
    --keep-matrix) KEEP_MATRIX=1; shift ;;
    -h|--help) usage ;;
    *) echo "unknown argument: $1" >&2; usage ;;
  esac
done

[[ -n "$CCRE_BED" && -n "$BIGWIG" && -n "$OUTFILE" ]] || usage
if [[ -z "$LABEL" ]]; then
  LABEL=$(basename "$CCRE_BED")
  LABEL="${LABEL%.*}"
fi

mkdir -p "$(dirname "$OUTFILE")"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# cCRE ids are suffixed "_<bin index>" and later split back apart on "_" to
# recover the id -- an id containing its own "_" would silently merge into
# the wrong group instead of erroring, so check for that up front.
bad_id=$(awk '{FS="\t"}{if ($4 ~ /_/) {print $4; exit}}' "$CCRE_BED")
if [[ -n "$bad_id" ]]; then
  echo "error: cCRE id '$bad_id' contains '_', which breaks id/bin-index reconstruction -- rename ids to remove '_'" >&2
  exit 1
fi

HALF=$((WINDOW / 2))

# ---------- 1. recenter each cCRE to a fixed-width window ----------
recentered="$TMP/recentered.bed4"
if [[ -n "$CHROM_SIZES" ]]; then
  awk -v half="$HALF" -v chromsizes="$CHROM_SIZES" '
    BEGIN{FS=OFS="\t"; while ((getline line < chromsizes) > 0) {split(line, a, "\t"); size[a[1]] = a[2]}}
    {
      center = int(($2 + $3) / 2)
      start = (center > half) ? center - half : 0
      end = center + half + 1
      if (($1 in size) && end > size[$1]) {end = size[$1]}
      print $1, start, end, $4
    }
  ' "$CCRE_BED" | sort -k1,1 -k2,2n > "$recentered"
else
  awk -v half="$HALF" '{FS=OFS="\t"}{
    center = int(($2 + $3) / 2)
    if (center > half) {print $1, center - half, center + half + 1, $4}
    else {print $1, 0, center + half + 1, $4}
  }' "$CCRE_BED" | sort -k1,1 -k2,2n > "$recentered"
fi
n_cres=$(wc -l < "$recentered")

# ---------- 2. split each window into N_BINS equal bins ----------
split_bed="$TMP/split.bed"
awk -v n="$N_BINS" '{FS=OFS="\t"}{
  l = $3 - $2
  l2 = l / n
  for (i = 1; i <= n; i++) {print $1, int($2 + l2 * (i - 1)), int($2 + l2 * i), $4"_"i}
}' "$recentered" > "$split_bed"

# ---------- 3. per-bin signal ----------
tab="$TMP/signal.tab"
bigWigAverageOverBed "$BIGWIG" "$split_bed" "$tab"

# ---------- 4. assemble the id x bin signal matrix ----------
matrix="$TMP/matrix.txt"
cut -f1,"$COLUMN" "$tab" \
  | awk '{FS=OFS="\t"}{split($1, a, "_"); print $1, a[1], a[2] + 0, $2}' \
  | sort -k2,2 -k3,3n \
  | awk -v n="$N_BINS" '{FS=OFS="\t"}{
      if ((NR - 1) % n == 0) {if (NR > 1) {printf "\n"}; printf $2}
      printf "\t"$4
    } END{printf "\n"}' > "$matrix"

if [[ "$KEEP_MATRIX" == "1" ]]; then
  cp "$matrix" "${OUTFILE}.matrix.txt"
fi

# ---------- 5. mean aggregation profile across all cCREs ----------
awk -v n="$N_BINS" -v label="$LABEL" '
  BEGIN{FS=OFS="\t"; for (i = 2; i <= n + 1; i++) {a[i] = 0}; nrows = 0}
  {nrows += 1; for (i = 2; i <= n + 1; i++) {a[i] += $i}}
  END{printf label; for (i = 2; i <= n + 1; i++) {printf "\t"a[i] / nrows}; printf "\n"}
' "$matrix" > "$OUTFILE"

echo "done: $n_cres cCREs aggregated into $N_BINS bins -> $OUTFILE"
