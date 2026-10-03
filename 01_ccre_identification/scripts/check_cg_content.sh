#!/bin/bash
set -euo pipefail

usage() {
  cat <<'USAGE' >&2
Usage: check_cg_content.sh -a CRE_BED4 -f GENOME_FASTA -o OUTFILE
                            [--monocg-bw BW] [--dicg-bw BW]

Check CG/GC sequence features of a given cCRE list: extracts each cCRE's
sequence (bedtools getfasta), computes its GC content and CpG
observed/expected ratio (CpG O/E = (CpG_count / length) /
((C_freq + G_freq) / 2)^2, a standard measure of CpG depletion/enrichment
relative to base composition), and optionally joins per-cCRE mono-CG and
di-CG signal from pre-built genome-wide CpG-density bigWig tracks.

Re-implementation (bash/awk/bedtools) of
0-7-4e_peak_sequence_feature.sh section 1 "CG feature", inlining
daily/calculate_CGdinucleotide.py and daily/count_GC_content.py as a single
awk pass (over properly concatenated multi-line sequences, same as the
GC-content script -- CpG O/E was computed line-by-line in the original,
which silently assumed one unwrapped line per record) so the script is
self-contained.

Requires bedtools on PATH (and bigWigAverageOverBed if --monocg-bw /
--dicg-bw are given). --monocg-bw/--dicg-bw are pre-built genome-wide
tracks (e.g. 0-7-4e's own prep step: scan the genome FASTA for mono-CG /
di-CG density, bedGraphToBigWig) -- building those is a separate,
genome-level prep step and out of scope here.

  -a CRE_BED4     cCRE BED, >=4 columns (chrom, start, end, id); id must be
                  unique
  -f GENOME_FASTA reference genome FASTA matching CRE_BED4's coordinates
  -o OUTFILE      output path for the merged per-cCRE CG feature table
  --monocg-bw BW  optional genome-wide mono-CG density bigWig; adds a
                  mono_cg column (mean signal over each cCRE)
  --dicg-bw BW    optional genome-wide di-CG density bigWig; adds a di_cg
                  column (mean signal over each cCRE)
USAGE
  exit 1
}

CRE_BED4=""
GENOME_FASTA=""
OUTFILE=""
MONOCG_BW=""
DICG_BW=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -a) CRE_BED4="$2"; shift 2 ;;
    -f) GENOME_FASTA="$2"; shift 2 ;;
    -o) OUTFILE="$2"; shift 2 ;;
    --monocg-bw) MONOCG_BW="$2"; shift 2 ;;
    --dicg-bw) DICG_BW="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "unknown argument: $1" >&2; usage ;;
  esac
done

[[ -n "$CRE_BED4" && -n "$GENOME_FASTA" && -n "$OUTFILE" ]] || usage

mkdir -p "$(dirname "$OUTFILE")"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# ---------- 1. extract cCRE sequences ----------
fasta="$TMP/cre.fa"
bedtools getfasta -nameOnly -fi "$GENOME_FASTA" -bed "$CRE_BED4" > "$fasta"

# ---------- 2. per-cCRE GC content + CpG observed/expected ----------
base_table="$TMP/base.txt"
awk '
  BEGIN{FS=OFS="\t"}
  /^>/ {
    if (id != "") {emit(id, seq)}
    id = substr($0, 2)
    sub(/\(.*\)$/, "", id)   # getfasta -nameOnly may still add "(+)"/"(-)" strand suffix
    seq = ""
    next
  }
  {seq = seq $0}
  END {if (id != "") {emit(id, seq)}}

  function emit(id, seq,   a, c, g, t, cg, len, total, real, expected, cpg_oe, gc_content) {
    seq = toupper(seq)
    len = length(seq)
    # gsub(pat, pat, seq) is a no-op replacement -- counts matches without
    # altering seq, so these calls can safely share one copy of seq
    a = gsub(/A/, "A", seq)
    c = gsub(/C/, "C", seq)
    g = gsub(/G/, "G", seq)
    t = gsub(/T/, "T", seq)
    cg = gsub(/CG/, "CG", seq)

    if (len > 0) {
      real = cg / len
      expected = ((c / len + g / len) / 2) ^ 2
    } else {
      expected = 0
    }
    cpg_oe = (expected == 0) ? "-1" : sprintf("%.2f", real / expected)

    total = a + c + g + t
    gc_content = (total > 0) ? sprintf("%.2f", (g + c) / total) : "0"

    print id, gc_content, cpg_oe, cg, c, g, len
  }
' "$fasta" > "$base_table"

cumulative="$base_table"

# ---------- 3. optional mono-CG / di-CG signal from pre-built bigWigs ----------
join_bigwig_signal() {
  # $1 = bigwig, $2 = column name -> appends a column to $cumulative
  #
  # NOTE: `colname` must be its own `local` statement before anything that
  # expands it -- in `local a=1 b=$a`, bash expands `$a` before the `local`
  # builtin actually creates `a`, so it would see any outer/unset `a`
  # instead of the new local value.
  local bw="$1" colname="$2"
  local tab="$TMP/${colname}.tab" new_cum="$TMP/cum.${colname}.txt"
  bigWigAverageOverBed "$bw" "$CRE_BED4" "$tab"
  awk -v tabfile="$tab" '
    BEGIN{FS=OFS="\t"}
    FNR==NR && FILENAME==tabfile {score[$1] = $5; next}
    {print $0, (($1 in score) ? score[$1] : "NA")}
  ' "$tab" "$cumulative" > "$new_cum"
  cumulative="$new_cum"
}

header="id\tgc_content\tcpg_oe\tnum_cg\tnum_c\tnum_g\tlength"
if [[ -n "$MONOCG_BW" ]]; then
  join_bigwig_signal "$MONOCG_BW" "mono_cg"
  header="${header}\tmono_cg"
fi
if [[ -n "$DICG_BW" ]]; then
  join_bigwig_signal "$DICG_BW" "di_cg"
  header="${header}\tdi_cg"
fi

echo -e "$header" > "$OUTFILE"
cat "$cumulative" >> "$OUTFILE"

echo "done: $(wc -l < "$base_table") cCREs -> $OUTFILE"
