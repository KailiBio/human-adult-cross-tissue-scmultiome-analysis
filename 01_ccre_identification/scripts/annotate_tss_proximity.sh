#!/bin/bash
# Group cCREs into TSS-overlap / TSS-proximal / TSS-distal by distance to
# the nearest GENCODE TSS.
#
# distance <= overlap_bp (default 200)   -> TSS-overlap
# distance <= proximal_bp (default 2000) -> TSS-proximal
# otherwise                              -> TSS-distal
#
# Output columns: all columns from -a (the cCRE bed), all columns from -b
# (the TSS bed), distance_to_tss, tss_category.
#
# Usage:
#   annotate_tss_proximity.sh -a final_ccres.bed -b gencode_tss.bed -o final_ccres.tss_annotated.tsv [-d 200] [-p 2000]
#
# Requires bedtools, sort, awk on PATH.

set -euo pipefail

usage() {
    echo "Usage: $0 -a <ccre.bed> -b <gencode_tss.bed> -o <output.tsv> [-d overlap_bp=200] [-p proximal_bp=2000]" >&2
    exit 1
}

OVERLAP_DIST=200
PROXIMAL_DIST=2000

while getopts "a:b:o:d:p:h" opt; do
    case "$opt" in
        a) CCRE_BED="$OPTARG" ;;
        b) TSS_BED="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        d) OVERLAP_DIST="$OPTARG" ;;
        p) PROXIMAL_DIST="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "${CCRE_BED:-}" || -z "${TSS_BED:-}" || -z "${OUTPUT:-}" ]]; then
    usage
fi

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

CCRE_SORTED="$TMPDIR/ccre.sorted.bed"
TSS_SORTED="$TMPDIR/tss.sorted.bed"

sort -k1,1 -k2,2n "$CCRE_BED" > "$CCRE_SORTED"
sort -k1,1 -k2,2n "$TSS_BED" > "$TSS_SORTED"

echo "computing distance to nearest GENCODE TSS" >&2

bedtools closest -a "$CCRE_SORTED" -b "$TSS_SORTED" -d -t first \
    | awk -v odist="$OVERLAP_DIST" -v pdist="$PROXIMAL_DIST" '
        BEGIN { FS = OFS = "\t" }
        {
            d = $NF
            if (d <= odist) { cat = "TSS-overlap" }
            else if (d <= pdist) { cat = "TSS-proximal" }
            else { cat = "TSS-distal" }
            print $0, cat
        }' > "$OUTPUT"

echo "TSS category counts:" >&2
awk -F'\t' '{print $NF}' "$OUTPUT" | sort | uniq -c
