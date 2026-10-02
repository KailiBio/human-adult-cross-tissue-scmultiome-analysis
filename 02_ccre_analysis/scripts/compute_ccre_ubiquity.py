#!/usr/bin/env python
"""Compute per-cCRE cross-tissue-celltype ubiquity.

Takes the final cCRE set and the per-cluster raw peaks it was built from
(both produced by 01_ccre_identification/scripts/call_ccres.py) and
produces, for each cCRE:

  1. the number of tissue-celltype clusters in which it is active
     (i.e. overlaps a raw peak called in that cluster)
  2. a binary cCRE x tissue-celltype matrix of that activity

Requires `bedtools` on PATH.
"""

import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

RAW_PEAK_COLS = [
    "chrom", "start", "end", "rawPeak", "score", "signal_value", "p_value", "q_value", "celltype",
]


def sort_bed(path: str, out_path: Path) -> None:
    with open(out_path, "w") as fh:
        subprocess.run(["sort", "-k1,1", "-k2,2n", path], stdout=fh, check=True)


def intersect_ccre_with_raw_peaks(ccre_bed: Path, raw_peaks_bed: Path, bedtools: str) -> pd.DataFrame:
    result = subprocess.run(
        [bedtools, "intersect", "-a", str(ccre_bed), "-b", str(raw_peaks_bed), "-wa", "-wb"],
        capture_output=True, text=True, check=True,
    )
    cols = ["chrom", "start", "end", "Peaks"] + ["raw_" + c for c in RAW_PEAK_COLS]
    rows = [line.split("\t") for line in result.stdout.splitlines() if line]
    return pd.DataFrame(rows, columns=cols)


def build_binary_matrix(overlaps: pd.DataFrame, all_peaks: list[str]) -> pd.DataFrame:
    active = overlaps[["Peaks", "raw_celltype"]].drop_duplicates()
    matrix = (
        active.assign(active=1)
        .pivot_table(index="Peaks", columns="raw_celltype", values="active", fill_value=0)
        .astype(np.int8)
        .reindex(all_peaks, fill_value=0)
    )
    matrix.index.name = "Peaks"
    return matrix


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--final-ccre-bed", required=True, help="final_ccres.bed from call_ccres.py (chrom, start, end, Peaks)")
    ap.add_argument("--raw-peaks-bed", required=True, help="raw_peaks.capped.bed from call_ccres.py")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--bedtools", default="bedtools")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_peaks = pd.read_csv(
        args.final_ccre_bed, sep="\t", header=None, names=["chrom", "start", "end", "Peaks"]
    )["Peaks"].tolist()

    ccre_sorted = outdir / "final_ccres.sorted.bed"
    raw_sorted = outdir / "raw_peaks.sorted.bed"
    sort_bed(args.final_ccre_bed, ccre_sorted)
    sort_bed(args.raw_peaks_bed, raw_sorted)

    print("intersecting final cCREs with per-cluster raw peaks")
    overlaps = intersect_ccre_with_raw_peaks(ccre_sorted, raw_sorted, args.bedtools)

    print("building binary cCRE x tissue-celltype matrix")
    matrix = build_binary_matrix(overlaps, all_peaks)
    matrix.to_csv(outdir / "ccre_celltype_binary_matrix.tsv", sep="\t")

    ubiquity = matrix.sum(axis=1).rename("n_active_celltypes").reset_index()
    ubiquity.to_csv(outdir / "ccre_ubiquity_counts.tsv", sep="\t", index=False)

    print(f"  {len(all_peaks)} cCREs x {matrix.shape[1]} tissue-celltype clusters")
    print(f"  median active celltypes per cCRE: {ubiquity['n_active_celltypes'].median()}")


if __name__ == "__main__":
    main()
