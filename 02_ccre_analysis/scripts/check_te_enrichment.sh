#!/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE' >&2
Usage: check_te_enrichment.sh -a CRE_INFO -t TE_BED -g CHROM_SIZES -o OUTDIR
                               [-i ID_COL] [-s GROUP_COL] [-e ENCODE_COL]
                               [-f FAMILY_COL] [-c CLASS_COL]
                               [--classes "DNA LINE LTR SINE Retroposon"]
                               [--sine-families "Alu MIR"]
                               [--ltr-families "ERVK ERV1 ERVL ERVL-MaLR"]
                               [--line-families "L1 L2"]

Check TE enrichment of cCREs: what fraction overlap a TE overall, by TE
class, and by TE family within SINE/LTR/LINE (+ an "other" bucket per
class), stratified by two grouping columns from CRE_INFO (default: a
category column and an ENCODE-overlap flag column) -- each compared against
that TE's background genome-wide coverage.

Each cCRE's midpoint (not its full interval) is used for overlap, so a wide
peak spanning several TE copies isn't double-counted.

Re-implementation (bash/awk/bedtools) of
0-7-4c_check_peak_specific_celltype.sh section 3 "enrichment with TE".
Requires `bedtools` on PATH.

  -a CRE_INFO      cCRE info table, tab-separated, no header (chrom, start,
                   end, ... with an id/group/encode-flag column each)
  -t TE_BED        RepeatMasker-derived TE BED, tab-separated, no header
  -g CHROM_SIZES   genome chrom.sizes (for background TE coverage)
  -o OUTDIR        output directory
  -i ID_COL        1-based column in CRE_INFO with the cCRE id (default: 4)
  -s GROUP_COL     1-based column in CRE_INFO with the first stratifying
                   label, e.g. tier/category (default: 7)
  -e ENCODE_COL    1-based column in CRE_INFO with the second stratifying
                   label, e.g. ENCODE-overlap flag (default: 9)
  -f FAMILY_COL    1-based column in TE_BED with the TE family, e.g.
                   'LINE/L1' (default: 7)
  -c CLASS_COL     1-based column in TE_BED with the TE class, e.g. 'LINE'
                   (default: 10)
  --classes        space-separated TE classes to test (default: "DNA LINE
                   LTR SINE Retroposon")
  --sine-families  space-separated SINE families to test individually;
                   everything else in SINE -> "other" (default: "Alu MIR")
  --ltr-families   as above for LTR (default: "ERVK ERV1 ERVL ERVL-MaLR")
  --line-families  as above for LINE (default: "L1 L2")
USAGE
  exit 1
}

ID_COL=4
GROUP_COL=7
ENCODE_COL=9
FAMILY_COL=7
CLASS_COL=10
CLASSES="DNA LINE LTR SINE Retroposon"
SINE_FAMILIES="Alu MIR"
LTR_FAMILIES="ERVK ERV1 ERVL ERVL-MaLR"
LINE_FAMILIES="L1 L2"
CRE_INFO=""
TE_BED=""
CHROM_SIZES=""
OUTDIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -a) CRE_INFO="$2"; shift 2 ;;
    -t) TE_BED="$2"; shift 2 ;;
    -g) CHROM_SIZES="$2"; shift 2 ;;
    -o) OUTDIR="$2"; shift 2 ;;
    -i) ID_COL="$2"; shift 2 ;;
    -s) GROUP_COL="$2"; shift 2 ;;
    -e) ENCODE_COL="$2"; shift 2 ;;
    -f) FAMILY_COL="$2"; shift 2 ;;
    -c) CLASS_COL="$2"; shift 2 ;;
    --classes) CLASSES="$2"; shift 2 ;;
    --sine-families) SINE_FAMILIES="$2"; shift 2 ;;
    --ltr-families) LTR_FAMILIES="$2"; shift 2 ;;
    --line-families) LINE_FAMILIES="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "unknown argument: $1" >&2; usage ;;
  esac
done

[[ -n "$CRE_INFO" && -n "$TE_BED" && -n "$CHROM_SIZES" && -n "$OUTDIR" ]] || usage

mkdir -p "$OUTDIR"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

BG_LOG="$OUTDIR/te_background_coverage.log"
: > "$BG_LOG"
CUM_COUNTER=0

genome_size=$(awk 'BEGIN{FS=OFS="\t";sum=0}{sum+=$2}END{print sum}' "$CHROM_SIZES")

log_background() {
  # $1 = label, $2 = bed subset -> appends label, bp, background pct to BG_LOG
  local label="$1" bed="$2" bp pct
  bp=$(awk 'BEGIN{FS=OFS="\t";sum=0}{sum+=($3-$2)}END{print sum+0}' "$bed")
  pct=$(awk -v bp="$bp" -v g="$genome_size" 'BEGIN{printf "%.6f", bp/g*100}')
  echo -e "${label}\t${bp}\t${pct}" | tee -a "$BG_LOG"
}

overlap_counts() {
  # $1 = TE subset bed -> group\tencode\tn_overlap on stdout
  bedtools intersect -a "$TMP/cre_center.bed" -b "$1" -wa \
    | cut -f4-6 | sort -u | cut -f2,3 | sort | uniq -c \
    | awk '{OFS="\t"}{print $2,$3,$1}'
}

append_counts() {
  # $1 = cumulative table so far, $2 = new overlap-count file (group, encode,
  # n_overlap); appends n_overlap,pct (of column 4, the "total" column) into
  # a fresh file and sets the global NEW_CUM to its path. Missing
  # (group,encode) rows in $2 are treated as 0 overlap.
  #
  # NOTE: must be called as a plain statement, not via `$(...)` -- command
  # substitution forks a subshell, so CUM_COUNTER's increment wouldn't
  # survive back to the caller and every call would silently reuse (and
  # truncate) the same output file.
  #
  # NOTE: the lookup-array file ($counts) can be completely empty (a class
  # or family with zero overlapping cCREs). Guard with FILENAME rather than
  # the usual `NR==FNR` -- if the first file has 0 records, FNR goes right
  # on equaling NR for every record of the *second* file too, which would
  # silently swallow the whole cumulative table into the lookup branch
  # instead of printing it.
  local cumulative="$1" counts="$2"
  CUM_COUNTER=$((CUM_COUNTER + 1))
  NEW_CUM="$TMP/cum_${CUM_COUNTER}.txt"
  awk -v countsfile="$counts" '{FS=OFS="\t"}{
    if(FNR==NR && FILENAME==countsfile){a[$1 SUBSEP $2]=$3; next}
    key=$1 SUBSEP $2
    n=(key in a)?a[key]:0
    print $0,n,n/$4*100
  }' "$counts" "$cumulative" > "$NEW_CUM"
}

te_class_bed() {
  # $1 = class, $2 = out path
  awk -v c="$CLASS_COL" -v cls="$1" '{FS=OFS="\t"}{if($c==cls){print $1,$2,$3}}' "$TE_BED" > "$2"
}

te_family_bed() {
  # $1 = full family name (e.g. 'SINE/Alu'), $2 = out path
  awk -v f="$FAMILY_COL" -v fam="$1" '{FS=OFS="\t"}{if($f==fam){print $1,$2,$3}}' "$TE_BED" > "$2"
}

te_other_bed() {
  # $1 = class, $2 = space-separated full family names to exclude, $3 = out path
  awk -v c="$CLASS_COL" -v f="$FAMILY_COL" -v cls="$1" -v excl="$2" '
    BEGIN{FS=OFS="\t"; n=split(excl,arr," "); for(i=1;i<=n;i++){skip[arr[i]]=1}}
    {if($c==cls && !($f in skip)){print $1,$2,$3}}
  ' "$TE_BED" > "$3"
}

# ---------- 0. cCRE midpoint bed: chrom, center, center+1, id, group, encode_flag ----------
awk -v id="$ID_COL" -v grp="$GROUP_COL" -v enc="$ENCODE_COL" \
  '{FS=OFS="\t"}{center=int(($2+$3)/2); print $1,center,center+1,$id,$grp,$enc}' \
  "$CRE_INFO" > "$TMP/cre_center.bed"

cut -f5,6 "$TMP/cre_center.bed" | sort | uniq -c \
  | awk '{OFS="\t"}{print $2,$3,$1}' > "$TMP/n_total.txt"

# ---------- 1. overall TE overlap ----------
echo "testing overlap with any TE"
overlap_counts "$TE_BED" > "$TMP/n_overlap_any.txt"
awk '{FS=OFS="\t"}{
  if(NR==FNR){a[$1 SUBSEP $2]=$3}
  else{print $0,a[$1 SUBSEP $2],$3/a[$1 SUBSEP $2]*100}
}' "$TMP/n_total.txt" "$TMP/n_overlap_any.txt" > "$OUTDIR/ccre_overlap_te.txt"
log_background "any_TE" "$TE_BED"

# ---------- 2. by TE class ----------
cumulative="$TMP/cum_0.txt"
cp "$OUTDIR/ccre_overlap_te.txt" "$cumulative"

for cls in $CLASSES; do
  echo "testing class $cls"
  cls_bed="$TMP/cls_${cls}.bed"
  te_class_bed "$cls" "$cls_bed"
  log_background "$cls" "$cls_bed"
  overlap_counts "$cls_bed" > "$TMP/counts.txt"
  append_counts "$cumulative" "$TMP/counts.txt"
  cumulative="$NEW_CUM"
done
cp "$cumulative" "$OUTDIR/ccre_overlap_te.by_class.txt"

# ---------- 3. by TE family within SINE / LTR / LINE (+ "other" per class) ----------
cumulative="$TMP/cum_0.txt"

process_family_group() {
  local cls="$1" fams="$2" fam full fam_bed other_bed full_names=""
  for fam in $fams; do
    full="${cls}/${fam}"
    full_names="${full_names} ${full}"
    echo "  testing family $full"
    fam_bed="$TMP/fam_${cls}_${fam}.bed"
    te_family_bed "$full" "$fam_bed"
    log_background "$full" "$fam_bed"
    overlap_counts "$fam_bed" > "$TMP/counts.txt"
    append_counts "$cumulative" "$TMP/counts.txt"
    cumulative="$NEW_CUM"
  done
  echo "  testing family ${cls}/other"
  other_bed="$TMP/fam_${cls}_other.bed"
  te_other_bed "$cls" "${full_names# }" "$other_bed"
  log_background "${cls}/other" "$other_bed"
  overlap_counts "$other_bed" > "$TMP/counts.txt"
  append_counts "$cumulative" "$TMP/counts.txt"
  cumulative="$NEW_CUM"
}

echo "testing class SINE families"
process_family_group "SINE" "$SINE_FAMILIES"
echo "testing class LTR families"
process_family_group "LTR" "$LTR_FAMILIES"
echo "testing class LINE families"
process_family_group "LINE" "$LINE_FAMILIES"

cp "$cumulative" "$OUTDIR/ccre_overlap_te.by_family.txt"

echo "done"
