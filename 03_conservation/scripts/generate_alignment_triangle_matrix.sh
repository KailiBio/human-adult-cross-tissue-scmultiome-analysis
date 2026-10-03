#!/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE' >&2
Usage: generate_alignment_triangle_matrix.sh -a CRE_BED4 -b ALIGN_BG_DIR -g CHROM_SIZES -o OUTDIR
                                              [-r RANDOM_BED4] [-i CRE_INFO] [-p PAR_BED]
                                              [--id-col N] [--group-col N] [--group2-col N]
                                              [--align-high 0.9] [--align-low 0.1] [--bw-dir DIR]

For each species with a per-species genome alignment bedGraph, computes
per-cCRE mean alignment score (via bigWigAverageOverBed) over the cCRE set
(and optionally a background/random region set), then tallies, per cCRE,
how many species have alignment >= --align-high ("align90") vs
<= --align-low ("align10"). These two counts are the x/y axes of the
"triangle plot" (named for the triangular point cloud align90+align10 <=
n_species produces) used to flag accelerated / non-conserved accessible
elements.

Re-implementation (bash/awk/bedtools/UCSC tools) of
0-8-1_generate_alignment_file.sh sections 0-3. Requires bedtools,
bedGraphToBigWig, and bigWigAverageOverBed on PATH.

  -a CRE_BED4     cCRE BED4 (chrom, start, end, id); id must be unique
  -b ALIGN_BG_DIR directory of per-species alignment bedGraphs, one
                  <species>.bg file per species
  -g CHROM_SIZES  genome chrom.sizes matching the alignment bedGraphs
  -o OUTDIR       output directory
  -r RANDOM_BED4  optional background/random-region BED4, scored the same
                  way for comparison (e.g. random genomic windows)
  -i CRE_INFO     optional cCRE info table (tab-separated, no header,
                  chrom/start/end in the first 3 columns so it can be
                  PAR-filtered) to stratify the cCRE triangle count table
                  by group/group2 (e.g. TSS-proximity tier, known/novel);
                  without it the count table is just (id, align90, align10)
  -p PAR_BED      optional pseudoautosomal-region BED; cCREs overlapping
                  it are excluded from the (stratified) triangle count
                  table
  --id-col        1-based column in CRE_INFO with the cCRE id (default: 4)
  --group-col     1-based column in CRE_INFO with the group label
                  (default: 5)
  --group2-col    1-based column in CRE_INFO with the group2 label
                  (default: 6)
  --align-high    alignment score threshold for "align90" (default: 0.9)
  --align-low     alignment score threshold for "align10" (default: 0.1)
  --bw-dir        directory to read/write per-species bigWigs (default:
                  ALIGN_BG_DIR/bw); existing <species>.bw files are reused
USAGE
  exit 1
}

ALIGN_HIGH=0.9
ALIGN_LOW=0.1
ID_COL=4
GROUP_COL=5
GROUP2_COL=6
BW_DIR=""
CRE_BED4=""
ALIGN_BG_DIR=""
CHROM_SIZES=""
OUTDIR=""
RANDOM_BED4=""
CRE_INFO=""
PAR_BED=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -a) CRE_BED4="$2"; shift 2 ;;
    -b) ALIGN_BG_DIR="$2"; shift 2 ;;
    -g) CHROM_SIZES="$2"; shift 2 ;;
    -o) OUTDIR="$2"; shift 2 ;;
    -r) RANDOM_BED4="$2"; shift 2 ;;
    -i) CRE_INFO="$2"; shift 2 ;;
    -p) PAR_BED="$2"; shift 2 ;;
    --id-col) ID_COL="$2"; shift 2 ;;
    --group-col) GROUP_COL="$2"; shift 2 ;;
    --group2-col) GROUP2_COL="$2"; shift 2 ;;
    --align-high) ALIGN_HIGH="$2"; shift 2 ;;
    --align-low) ALIGN_LOW="$2"; shift 2 ;;
    --bw-dir) BW_DIR="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "unknown argument: $1" >&2; usage ;;
  esac
done

[[ -n "$CRE_BED4" && -n "$ALIGN_BG_DIR" && -n "$CHROM_SIZES" && -n "$OUTDIR" ]] || usage
[[ -z "$BW_DIR" ]] && BW_DIR="$ALIGN_BG_DIR/bw"

mkdir -p "$OUTDIR" "$BW_DIR" "$OUTDIR/per_species"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# ---------- 0. species list ----------
species_list="$OUTDIR/species_list.txt"
: > "$species_list"
for bg in "$ALIGN_BG_DIR"/*.bg; do
  species=$(basename "$bg" .bg)
  echo "$species" >> "$species_list"
done
n_species=$(wc -l < "$species_list")
echo "found $n_species species in $ALIGN_BG_DIR"

# ---------- 1. per-species alignment bigWig (reused if already built) ----------
while read -r species; do
  bw="$BW_DIR/${species}.bw"
  if [[ -f "$bw" ]]; then
    echo "reusing $bw"
  else
    echo "building $bw"
    bedGraphToBigWig "$ALIGN_BG_DIR/${species}.bg" "$CHROM_SIZES" "$bw"
  fi
done < "$species_list"

# ---------- 2. per-species alignment score over the cCRE set (and background, if given) ----------
CUM_COUNTER=0

score_one_set() {
  # $1 = bed4, $2 = label (e.g. "cre"/"background") -> writes $OUTDIR/per_species/<species>_alignment.<label>.tab per species
  local bed4="$1" label="$2" species bw
  while read -r species; do
    echo "scoring $label: $species"
    bw="$BW_DIR/${species}.bw"
    bigWigAverageOverBed "$bw" "$bed4" "$OUTDIR/per_species/${species}_alignment.${label}.tab"
  done < "$species_list"
}

build_score_matrix() {
  # $1 = bed4, $2 = label -> writes $OUTDIR/alignment_score_matrix.<label>.txt (id, one column per species)
  #
  # NOTE: `label` must be its own `local` statement before anything that
  # expands it -- in `local a=1 b=$a`, bash expands `$a` before the `local`
  # builtin actually creates `a`, so it would see any outer/unset `a`
  # instead of the new local value.
  local bed4="$1" label="$2"
  local out="$OUTDIR/alignment_score_matrix.${label}.txt"
  local header="id" tmp_paste="$TMP/paste.${label}.txt" species
  : > "$tmp_paste"
  while read -r species; do
    header="${header}\t${species}"
    awk '{print $5}' "$OUTDIR/per_species/${species}_alignment.${label}.tab" > "$TMP/col.${label}.${species}.txt"
  done < "$species_list"

  cut -f4 "$bed4" > "$TMP/ids.${label}.txt"
  cols=("$TMP/ids.${label}.txt")
  while read -r species; do
    cols+=("$TMP/col.${label}.${species}.txt")
  done < "$species_list"

  echo -e "$header" > "$out"
  paste "${cols[@]}" >> "$out"
}

# NOTE: append_triangle_counts must be called as a plain statement, not via
# `$(...)` -- command substitution forks a subshell, so CUM_COUNTER's
# increment wouldn't survive back to the caller.
#
# NOTE: guards the lookup-array build with FILENAME rather than the usual
# `NR==FNR` -- a species whose alignment tab has 0 matching rows (e.g. a
# cCRE set outside that species' alignment coverage) would otherwise make
# FNR==NR true for every row of the cumulative table too, silently
# swallowing it instead of passing it through.
append_triangle_counts() {
  # $1 = cumulative table so far, $2 = species tab (name, size, covered, sum,
  # mean0, mean), $3 = align90 column (1-based), $4 = align10 column
  local cumulative="$1" tab="$2" c90="$3" c10="$4"
  CUM_COUNTER=$((CUM_COUNTER + 1))
  NEW_CUM="$TMP/triangle_cum_${CUM_COUNTER}.txt"
  awk -v tabfile="$tab" -v hi="$ALIGN_HIGH" -v lo="$ALIGN_LOW" -v c90="$c90" -v c10="$c10" '
    BEGIN{FS=OFS="\t"}
    FNR==NR && FILENAME==tabfile {score[$1]=$5; next}
    {
      if ($1 in score) {
        s = score[$1] + 0
        if (s >= hi) {$c90 = $c90 + 1}
        else if (s <= lo) {$c10 = $c10 + 1}
      }
      print $0
    }
  ' "$tab" "$cumulative" > "$NEW_CUM"
}

build_triangle_counts_cre() {
  local cumulative species
  if [[ -n "$CRE_INFO" ]]; then
    local info="$CRE_INFO"
    if [[ -n "$PAR_BED" ]]; then
      info="$TMP/cre_info.rmPAR.txt"
      bedtools intersect -a "$CRE_INFO" -b "$PAR_BED" -wa -v > "$info"
    fi
    echo -e "id\tgroup\tgroup2\talign90\talign10" > "$TMP/triangle_cum_0.txt"
    awk -v id="$ID_COL" -v g="$GROUP_COL" -v g2="$GROUP2_COL" \
      '{FS=OFS="\t"}{print $id,$g,$g2,0,0}' "$info" >> "$TMP/triangle_cum_0.txt"
    cumulative="$TMP/triangle_cum_0.txt"
    c90=4; c10=5
  else
    echo -e "id\talign90\talign10" > "$TMP/triangle_cum_0.txt"
    awk '{FS=OFS="\t"}{print $4,0,0}' "$CRE_BED4" >> "$TMP/triangle_cum_0.txt"
    cumulative="$TMP/triangle_cum_0.txt"
    c90=2; c10=3
  fi

  while read -r species; do
    echo "tallying align90/align10: $species"
    append_triangle_counts "$cumulative" "$OUTDIR/per_species/${species}_alignment.cre.tab" "$c90" "$c10"
    cumulative="$NEW_CUM"
  done < "$species_list"

  cp "$cumulative" "$OUTDIR/alignment_triangle_counts.cre.txt"
}

build_triangle_counts_background() {
  local cumulative species
  echo -e "id\talign90\talign10" > "$TMP/triangle_bg_0.txt"
  awk '{FS=OFS="\t"}{print $4,0,0}' "$RANDOM_BED4" >> "$TMP/triangle_bg_0.txt"
  cumulative="$TMP/triangle_bg_0.txt"

  while read -r species; do
    echo "tallying align90/align10 (background): $species"
    append_triangle_counts "$cumulative" "$OUTDIR/per_species/${species}_alignment.background.tab" 2 3
    cumulative="$NEW_CUM"
  done < "$species_list"

  cp "$cumulative" "$OUTDIR/alignment_triangle_counts.background.txt"
}

score_one_set "$CRE_BED4" "cre"
build_score_matrix "$CRE_BED4" "cre"
build_triangle_counts_cre

if [[ -n "$RANDOM_BED4" ]]; then
  score_one_set "$RANDOM_BED4" "background"
  build_score_matrix "$RANDOM_BED4" "background"
  build_triangle_counts_background
fi

echo "done"
