#!/usr/bin/env python
"""Call a final consensus cCRE set from per-cluster scATAC fragments.

Pipeline:
  1. Call raw MACS3 peaks per tissue-celltype cluster (snapatac2).
  2. Keep only clusters with enough sequencing depth (cells or fragments),
     cap each kept cluster to its top-N peaks by q-value, and merge into a
     consensus peak set.
  3. Intersect consensus peaks against the raw per-cluster peaks and keep
     the best (max) q-value supporting each consensus peak, then threshold
     on q-value to produce the final cCRE set.

Requires `bedtools` on PATH.
"""

import argparse
import subprocess
from pathlib import Path

import pandas as pd
import scanpy as sc
import snapatac2 as snap


def read_chrom_sizes(path: str) -> dict[str, int]:
    chrom_sizes = {}
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            chrom, size = parts
            chrom_sizes[chrom] = int(size)
    return chrom_sizes


def call_raw_peaks(adata, groupby: str, blacklist: str, n_jobs: int) -> None:
    """Step 1: call raw MACS3 peaks per cluster, stored in adata.uns['macs3']."""
    snap.tl.macs3(adata, groupby=groupby, blacklist=blacklist, n_jobs=n_jobs)


def compute_cluster_depth(adata, groupby: str) -> pd.DataFrame:
    """Per-cluster cell count, total fragments, and raw peak count."""
    depth = adata.obs.groupby(groupby).agg(
        cell_count=("n_fragment", "count"),
        n_fragment=("n_fragment", "sum"),
    )
    peak_counts = pd.Series(
        {k: len(df) for k, df in adata.uns["macs3"].items()}, name="peak_count"
    )
    peak_counts.index.name = groupby
    return depth.join(peak_counts)


def select_deep_clusters(
    depth: pd.DataFrame, min_cells: int, min_fragments: int
) -> list[str]:
    """Clusters with enough sequencing depth to call reliable peaks from."""
    mask = (depth["cell_count"] >= min_cells) | (depth["n_fragment"] >= min_fragments)
    return depth.index[mask].tolist()


def cap_and_merge_peaks(
    macs3: dict[str, pd.DataFrame],
    clusters: list[str],
    top_n_peaks: int,
    chrom_sizes: dict[str, int],
    merge_half_width: int,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """Step 2: cap each deep-enough cluster to its top-N peaks by q-value,
    then merge across clusters into a consensus peak set."""
    capped = {
        cluster: macs3[cluster].nlargest(top_n_peaks, "q_value")
        for cluster in clusters
        if cluster in macs3
    }
    merged = snap.tl.merge_peaks(capped, chrom_sizes, half_width=merge_half_width)
    consensus = pd.DataFrame(merged)
    consensus.columns = merged.columns  # pd.DataFrame() on the polars result drops column names
    return capped, consensus


def write_raw_peaks_bed(capped: dict[str, pd.DataFrame], path: Path) -> None:
    cols = ["chrom", "start", "end", "score", "signal_value", "p_value", "q_value"]
    rows = []
    for cluster, df in capped.items():
        if df is None or len(df) == 0:
            continue
        sub = df[cols].copy()
        sub["rawPeak"] = sub["chrom"] + ":" + sub["start"].astype(str) + "-" + sub["end"].astype(str)
        sub["celltype"] = cluster
        rows.append(sub)
    all_peaks = pd.concat(rows, ignore_index=True)
    all_peaks[["chrom", "start", "end", "rawPeak", "score", "signal_value", "p_value", "q_value", "celltype"]].to_csv(
        path, sep="\t", header=False, index=False
    )


def write_consensus_bed(consensus: pd.DataFrame, path: Path) -> pd.DataFrame:
    region = pd.DataFrame(consensus["Peaks"])
    region[["chrom", "pos"]] = region["Peaks"].str.split(":", expand=True)
    region[["start", "end"]] = region["pos"].str.split("-", expand=True)
    region["start"] = region["start"].astype(int)
    region["end"] = region["end"].astype(int)
    region[["chrom", "start", "end", "Peaks"]].to_csv(path, sep="\t", header=False, index=False)
    return region


def intersect_with_bedtools(consensus_bed: Path, raw_bed: Path, out: Path, bedtools: str) -> None:
    with open(out, "w") as fh:
        subprocess.run(
            [bedtools, "intersect", "-a", str(consensus_bed), "-b", str(raw_bed), "-wa", "-wb"],
            stdout=fh,
            check=True,
        )


def pick_best_qvalue(overlaps_path: Path) -> pd.DataFrame:
    """Step 3: for each consensus peak, keep the raw-peak overlap with the
    highest q-value (MACS3 reports -log10(q), so larger = more significant)."""
    cols = [
        "chrom", "start", "end", "Peaks",
        "b_chrom", "b_start", "b_end",
        "rawPeak", "score", "signal_value", "p_value", "q_value", "celltype",
    ]
    overlaps = pd.read_csv(overlaps_path, sep="\t", header=None, names=cols)
    overlaps = overlaps[overlaps["q_value"].notnull()]
    best_idx = overlaps.groupby("Peaks")["q_value"].idxmax()
    best = overlaps.loc[
        best_idx,
        ["chrom", "start", "end", "Peaks", "rawPeak", "score", "signal_value", "p_value", "q_value", "celltype"],
    ]
    return best.reset_index(drop=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-h5ad", required=True, help="Fragment-level AnnData (obs needs --groupby column and 'n_fragment')")
    ap.add_argument("--groupby", default="celltype_glue", help="obs column defining tissue-celltype clusters")
    ap.add_argument("--blacklist", required=True, help="BED file of blacklisted regions (e.g. ENCODE hg38 blacklist)")
    ap.add_argument("--chrom-sizes", required=True, help="Two-column chrom.sizes file")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-cells", type=int, default=300, help="Min cells per cluster to call peaks from")
    ap.add_argument("--min-fragments", type=int, default=1_000_000, help="Min total fragments per cluster to call peaks from")
    ap.add_argument("--top-n-peaks", type=int, default=180_000, help="Max peaks kept per cluster, ranked by q-value, before merging")
    ap.add_argument("--merge-half-width", type=int, default=150, help="Half-width (bp) used when merging per-cluster peaks into consensus regions")
    ap.add_argument("--qvalue-threshold", type=float, default=4.0, help="Keep consensus peaks whose best supporting q-value (-log10 q) exceeds this")
    ap.add_argument("--n-jobs", type=int, default=24)
    ap.add_argument("--bedtools", default="bedtools")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    chrom_sizes = read_chrom_sizes(args.chrom_sizes)

    print(f"[1/3] loading {args.input_h5ad} and calling raw peaks per {args.groupby}")
    adata = sc.read_h5ad(args.input_h5ad)
    call_raw_peaks(adata, groupby=args.groupby, blacklist=args.blacklist, n_jobs=args.n_jobs)
    adata.write(outdir / "atac_raw_peaks.h5ad")

    print("[2/3] checking cluster depth and merging consensus peaks")
    depth = compute_cluster_depth(adata, groupby=args.groupby)
    depth.to_csv(outdir / "cluster_depth_and_peak_counts.tsv", sep="\t")

    deep_clusters = select_deep_clusters(depth, args.min_cells, args.min_fragments)
    print(f"  {len(deep_clusters)}/{len(depth)} clusters pass depth filter "
          f"(cells>={args.min_cells} or fragments>={args.min_fragments})")

    capped, consensus = cap_and_merge_peaks(
        adata.uns["macs3"], deep_clusters, args.top_n_peaks, chrom_sizes, args.merge_half_width
    )
    print(f"  {consensus.shape[0]} consensus peaks before q-value filtering")

    consensus_bed = outdir / "consensus_peaks.bed"
    raw_bed = outdir / "raw_peaks.capped.bed"
    overlaps_bed = outdir / "consensus_vs_raw_overlaps.bed"

    write_consensus_bed(consensus, consensus_bed)
    write_raw_peaks_bed(capped, raw_bed)

    print("[3/3] intersecting with bedtools and filtering by best q-value")
    intersect_with_bedtools(consensus_bed, raw_bed, overlaps_bed, args.bedtools)
    best = pick_best_qvalue(overlaps_bed)

    final = best[best["q_value"] > args.qvalue_threshold].copy()
    best.to_csv(outdir / "consensus_peaks.best_qvalue.tsv", sep="\t", index=False)
    final[["chrom", "start", "end", "Peaks"]].sort_values(["chrom", "start"]).to_csv(
        outdir / "final_ccres.bed", sep="\t", header=False, index=False
    )
    final.to_csv(outdir / "final_ccres.tsv", sep="\t", index=False)

    print(f"  {len(best)} consensus peaks scored, {len(final)} pass "
          f"q-value > {args.qvalue_threshold} -> final_ccres.bed")


if __name__ == "__main__":
    main()
