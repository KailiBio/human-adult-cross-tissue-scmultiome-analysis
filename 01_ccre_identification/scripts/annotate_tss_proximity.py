#!/usr/bin/env python
"""Group cCREs into TSS-overlap / TSS-proximal / TSS-distal by distance to GENCODE TSS.

For each cCRE, finds the distance to the nearest GENCODE TSS (0 if the cCRE
overlaps a TSS) via `bedtools closest -d`, then buckets cCREs by that
distance:

  distance <= --overlap-distance   -> TSS-overlap   (default: 200 bp)
  distance <= --proximal-distance  -> TSS-proximal   (default: 2,000 bp)
  otherwise                        -> TSS-distal

Requires `bedtools` on PATH.
"""

import argparse
import subprocess
import tempfile
from pathlib import Path

import pandas as pd


def sort_bed(path: str, out_path: Path) -> None:
    with open(out_path, "w") as fh:
        subprocess.run(["sort", "-k1,1", "-k2,2n", path], stdout=fh, check=True)


def count_columns(path: Path) -> int:
    with open(path) as fh:
        return len(fh.readline().rstrip("\n").split("\t"))


def closest_tss_distance(ccre_bed: Path, tss_bed: Path, bedtools: str) -> pd.DataFrame:
    """One row per cCRE: its coordinates/id plus distance to the nearest TSS."""
    n_ccre_cols = count_columns(ccre_bed)
    n_tss_cols = count_columns(tss_bed)

    result = subprocess.run(
        [bedtools, "closest", "-a", str(ccre_bed), "-b", str(tss_bed), "-d", "-t", "first"],
        capture_output=True, text=True, check=True,
    )
    rows = [line.split("\t") for line in result.stdout.splitlines() if line]

    ccre_cols = [f"col{i}" for i in range(n_ccre_cols)]
    ccre_cols[:4] = ["chrom", "start", "end", "Peaks"][: min(4, n_ccre_cols)]
    tss_cols = [f"tss_col{i}" for i in range(n_tss_cols)]
    tss_cols[:3] = ["tss_chrom", "tss_start", "tss_end"][: min(3, n_tss_cols)]
    cols = ccre_cols + tss_cols + ["distance_to_tss"]

    df = pd.DataFrame(rows, columns=cols)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["distance_to_tss"] = df["distance_to_tss"].astype(int)
    return df


def classify_tss_proximity(
    df: pd.DataFrame, overlap_distance: int, proximal_distance: int
) -> pd.DataFrame:
    def bucket(d: int) -> str:
        if d <= overlap_distance:
            return "TSS-overlap"
        if d <= proximal_distance:
            return "TSS-proximal"
        return "TSS-distal"

    df = df.copy()
    df["tss_category"] = df["distance_to_tss"].apply(bucket)
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ccre-bed", required=True, help="Final cCRE BED (chrom, start, end, id, ...)")
    ap.add_argument("--tss-bed", required=True, help="GENCODE TSS positions BED (chrom, start, end, ...)")
    ap.add_argument("--output", required=True, help="Output TSV path with distance + category columns")
    ap.add_argument("--overlap-distance", type=int, default=200, help="Max distance (bp) to call TSS-overlap")
    ap.add_argument("--proximal-distance", type=int, default=2000, help="Max distance (bp) to call TSS-proximal")
    ap.add_argument("--bedtools", default="bedtools")
    args = ap.parse_args()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        ccre_sorted = Path(tmpdir) / "ccre.sorted.bed"
        tss_sorted = Path(tmpdir) / "tss.sorted.bed"
        sort_bed(args.ccre_bed, ccre_sorted)
        sort_bed(args.tss_bed, tss_sorted)

        print("computing distance to nearest GENCODE TSS")
        distances = closest_tss_distance(ccre_sorted, tss_sorted, args.bedtools)

    annotated = classify_tss_proximity(distances, args.overlap_distance, args.proximal_distance)
    annotated.to_csv(output, sep="\t", index=False)

    counts = annotated["tss_category"].value_counts()
    print(counts.to_string())


if __name__ == "__main__":
    main()
