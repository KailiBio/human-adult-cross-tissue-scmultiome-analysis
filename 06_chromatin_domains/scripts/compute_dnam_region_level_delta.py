#!/usr/bin/env python
"""Per-ADR-region tumor-vs-normal DNA methylation delta, aggregated into
subgroups by equal-region weighting (not equal-probe weighting).

Takes per-sample, per-ADR-region mean beta (long format) and a
region_id -> region_type/subgroup map (an ADR BED), computes
delta = mean(Tumor) - mean(Normal) independently for each individual ADR
region within each TCGA project, then aggregates those per-region deltas
within each subgroup (mean or median across regions) -- unlike pooling all
probes in a subgroup together (compute_dnam_tumor_normal_delta.py), this
gives every ADR region equal weight regardless of how many probes it
contains.

The per-sample per-region mean beta input is a separate, heavier data-prep
step and out of scope here (raw IDATs -> beta values -> per-sample,
per-region mean).

Consolidated from methyl_work/scripts/analyze_region_level.py
(1-3-5b_check_TCGA_DNAme.ipynb section 3).
"""

import argparse
import glob
from pathlib import Path

import pandas as pd


def region_deltas(reg: pd.DataFrame, value_col: str, region_map: pd.DataFrame, min_samples: int) -> pd.DataFrame:
    """project x region_id delta (mean Tumor - mean Normal), per region."""
    g = reg.dropna(subset=[value_col]).groupby(["project", "region_id", "tissue_type"])[value_col]
    stats = g.agg(["mean", "size"]).unstack("tissue_type")
    stats.columns = [f"{a}_{b}" for a, b in stats.columns]

    out = stats.reset_index()
    keep = (out.get("size_Tumor", 0) >= min_samples) & (out.get("size_Normal", 0) >= min_samples)
    out = out[keep].copy()
    out = out.rename(columns={
        "mean_Tumor": "mean_tumor", "mean_Normal": "mean_normal",
        "size_Tumor": "n_tumor", "size_Normal": "n_normal",
    })
    out["delta_tumor_minus_normal"] = out["mean_tumor"] - out["mean_normal"]
    out = out[["project", "region_id", "n_tumor", "n_normal", "mean_tumor", "mean_normal", "delta_tumor_minus_normal"]]
    out = out.merge(region_map, on="region_id", how="left")
    return out


def subgroup_agg_of_region_deltas(region_delta_df: pd.DataFrame, how: str, order) -> pd.DataFrame:
    piv = region_delta_df.groupby(["project", "subgroup"])["delta_tumor_minus_normal"].agg(how).unstack()
    order = [s for s in order if s in piv.columns] if order else sorted(piv.columns)
    return piv[order]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--adr-bed",
        required=True,
        help="ADR BED, tab-separated, no header: chr, start, end, region_id, celltype_desc, numeric_code, region_type, subgroup",
    )
    ap.add_argument(
        "--region-detail-glob",
        required=True,
        action="append",
        help="glob pattern matching per-sample per-region parquet files (repeatable); columns: sample_uuid, region_id, mean_beta, n",
    )
    ap.add_argument(
        "--sample-summary",
        required=True,
        help="per-sample summary CSV with sample_uuid, project, tissue_type, global_mean (e.g. all_projects_summary.csv from compute_dnam_tumor_normal_delta.py)",
    )
    ap.add_argument("--min-samples", type=int, default=2, help="minimum tumor/normal samples required to test a region in a project (default: 2)")
    ap.add_argument("--subgroup-order", default="", help="comma-separated preferred column order for the subgroup heatmap tables (default: discovery order)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    adr = pd.read_csv(
        args.adr_bed, sep="\t", header=None,
        names=["chr", "start", "end", "region_id", "celltype_desc", "numeric_code", "region_type", "subgroup"],
    )
    region_map = adr[["region_id", "region_type", "subgroup"]].drop_duplicates()

    region_files = [f for pattern in args.region_detail_glob for f in glob.glob(pattern)]
    if not region_files:
        raise SystemExit("no files matched --region-detail-glob")
    reg = pd.concat([pd.read_parquet(f) for f in region_files], ignore_index=True)
    reg = reg.drop_duplicates(subset=["sample_uuid", "region_id"])
    print(f"per-sample per-region rows loaded: {reg.shape[0]}")

    meta = pd.read_csv(args.sample_summary, usecols=["sample_uuid", "project", "tissue_type", "global_mean"])
    meta = meta.drop_duplicates(subset="sample_uuid")

    reg = reg.merge(meta, on="sample_uuid", how="inner")
    reg = reg.merge(region_map, on="region_id", how="inner")
    reg["mean_beta_globalnorm"] = reg["mean_beta"] - reg["global_mean"]
    print(f"per-sample per-region rows after merge: {reg.shape[0]} ({reg['sample_uuid'].nunique()} samples)")

    sg_order = [s for s in args.subgroup_order.split(",") if s]

    region_delta_raw = region_deltas(reg, "mean_beta", region_map, args.min_samples)
    region_delta_raw.to_csv(outdir / "region_level_delta.csv", index=False)
    print(f"region_level_delta.csv: {region_delta_raw.shape[0]} rows")

    region_delta_gn = region_deltas(reg, "mean_beta_globalnorm", region_map, args.min_samples)
    region_delta_gn.to_csv(outdir / "region_level_delta_globalnorm.csv", index=False)
    print(f"region_level_delta_globalnorm.csv: {region_delta_gn.shape[0]} rows")

    subgroup_agg_of_region_deltas(region_delta_raw, "mean", sg_order).to_csv(outdir / "heatmap_table_subgroup_delta_by_region.csv")
    subgroup_agg_of_region_deltas(region_delta_gn, "mean", sg_order).to_csv(outdir / "heatmap_table_subgroup_delta_by_region_globalnorm.csv")
    subgroup_agg_of_region_deltas(region_delta_raw, "median", sg_order).to_csv(outdir / "heatmap_table_subgroup_delta_by_region_median.csv")
    subgroup_agg_of_region_deltas(region_delta_gn, "median", sg_order).to_csv(outdir / "heatmap_table_subgroup_delta_by_region_median_globalnorm.csv")

    print("done")


if __name__ == "__main__":
    main()
