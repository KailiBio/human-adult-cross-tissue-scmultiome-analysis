#!/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE' >&2
Usage: get_chromatin_domain_bins.sh -g CHROM_SIZES -m MAPPABILITY_BW -o OUTFILE
                                     [-w BIN_SIZE] [-t THRESHOLD]
                                     [-b BLACKLIST_BED ...]

Tile the genome into fixed-size bins (default 100kb) as the base unit for
calling chromatin domains, then filter out bins with low mappability
and/or overlapping a blacklist, so downstream domain-calling only
considers "clean" genomic bins.

Re-implementation (bash/awk/bedtools) of 1-3-1d_ATACdepletedRegions.sh
section 0 "get clean 100kbs". Requires bedtools on PATH, and
bigWigAverageOverBed for mappability scoring.

  -g CHROM_SIZES    genome chrom.sizes
  -m MAPPABILITY_BW genome-wide mappability bigWig (e.g. k50 Umap
                    multi-track mappability)
  -o OUTFILE        output path for the clean bin BED4 (chrom, start, end,
                    "chrom:start-end")
  -w BIN_SIZE       bin width in bp (default: 100000)
  -t THRESHOLD      minimum mean mappability to keep a bin (default: 0.5)
  -b BLACKLIST_BED  a blacklist-like BED to exclude bins overlapping it
                    (e.g. ENCODE blacklist, unusual regions, GRC
                    exclusions); repeatable -- without any, no blacklist
                    filtering is applied
USAGE
  exit 1
}

BIN_SIZE=100000
THRESHOLD=0.5
CHROM_SIZES=""
MAPPABILITY_BW=""
OUTFILE=""
BLACKLIST_BEDS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -g) CHROM_SIZES="$2"; shift 2 ;;
    -m) MAPPABILITY_BW="$2"; shift 2 ;;
    -o) OUTFILE="$2"; shift 2 ;;
    -w) BIN_SIZE="$2"; shift 2 ;;
    -t) THRESHOLD="$2"; shift 2 ;;
    -b) BLACKLIST_BEDS+=("$2"); shift 2 ;;
    -h|--help) usage ;;
    *) echo "unknown argument: $1" >&2; usage ;;
  esac
done

[[ -n "$CHROM_SIZES" && -n "$MAPPABILITY_BW" && -n "$OUTFILE" ]] || usage

mkdir -p "$(dirname "$OUTFILE")"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# ---------- 1. tile the genome into bins ----------
bins="$TMP/genome_bins.bed4"
bedtools makewindows -g "$CHROM_SIZES" -w "$BIN_SIZE" \
  | awk '{FS=OFS="\t"}{print $1, $2, $3, $1":"$2"-"$3}' > "$bins"
n_total=$(wc -l < "$bins")

# ---------- 2. filter by mappability ----------
mappability_tab="$TMP/mappability.tab"
bigWigAverageOverBed "$MAPPABILITY_BW" "$bins" "$mappability_tab"

mappable="$TMP/mappable.bed4"
awk -v thr="$THRESHOLD" '{FS=OFS="\t"}{
  if (NR == FNR) {if ($5 > thr) {keep[$1] = 1}}
  else {if ($4 in keep) {print $0}}
}' "$mappability_tab" "$bins" > "$mappable"
n_mappable=$(wc -l < "$mappable")

# ---------- 3. remove blacklisted bins ----------
current="$mappable"
if [[ ${#BLACKLIST_BEDS[@]} -gt 0 ]]; then
  combined_blacklist="$TMP/blacklist.combined.bed"
  cat "${BLACKLIST_BEDS[@]}" | cut -f1-3 | sort -u | sort -k1,1 -k2,2n > "$combined_blacklist"
  clean="$TMP/clean.bed4"
  bedtools intersect -a "$current" -b "$combined_blacklist" -wa -v > "$clean"
  current="$clean"
fi

cp "$current" "$OUTFILE"
n_clean=$(wc -l < "$OUTFILE")

echo "tiled genome into $n_total ${BIN_SIZE}bp bins"
echo "$n_mappable passed mappability > $THRESHOLD"
echo "$n_clean remain after blacklist filtering -> $OUTFILE"
