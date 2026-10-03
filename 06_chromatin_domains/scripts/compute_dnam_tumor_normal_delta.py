#!/usr/bin/env python
"""Tumor-vs-normal DNA methylation delta per TCGA project, pooled across
probes in each ADR region-type / fine-grained subgroup.

Combines per-batch TCGA methylation summaries (one row per sample; a
region-type/subgroup column holds that sample's probe-pooled mean beta
across all probes in that category) into Mann-Whitney U tumor-vs-normal
tests and delta = mean(Tumor) - mean(Normal), per TCGA project. Computed
both on the raw region/subgroup mean beta, and background-normalized
(each sample's value minus its own whole-array global mean beta, to
control for bulk/global hypomethylation that varies by cancer type).

Input batch summaries are the per-sample output of a methylation
array-processing step (raw IDATs -> beta values -> joined against a
probe-to-ADR-region map -> per-sample mean beta per region type/subgroup;
a separate, heavier data-prep step and out of scope here): columns
sample_uuid, project, tissue_type, global_mean, `{region_type}_mean`/`_n`
per region type, `sg_{subgroup}_mean`/`_n` per subgroup.

Consolidated from methyl_work/scripts/analyze_full.py and
analyze_global.py (1-3-5b_check_TCGA_DNAme.ipynb sections 1-2).
"""

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def compute_stats(df: pd.DataFrame, value_cols, col_fn, label_col_name: str) -> pd.DataFrame:
    """Per-project Mann-Whitney U + delta for each label in value_cols.

    col_fn(label) returns the column name holding that label's per-sample
    mean beta (raw or already background-normalized).
    """
    records = []
    for project, g in df.groupby("project"):
        tumor = g[g.tissue_type == "Tumor"]
        normal = g[g.tissue_type == "Normal"]
        for label in value_cols:
            col = col_fn(label)
            if col not in df.columns:
                continue
            t = tumor[col].dropna()
            n = normal[col].dropna()
            if len(t) < 2 or len(n) < 2:
                continue
            _, p = mannwhitneyu(t, n, alternative="two-sided")
            records.append({
                "project": project, label_col_name: label,
                "n_tumor": len(t), "n_normal": len(n),
                "mean_tumor": t.mean(), "mean_normal": n.mean(),
                "delta_tumor_minus_normal": t.mean() - n.mean(),
                "pval": p,
            })
    out = pd.DataFrame(records)
    if not out.empty:
        out["neglog10p"] = -np.log10(out["pval"].clip(lower=1e-300))
    return out


def heatmap_table(stats: pd.DataFrame, label_col_name: str, order) -> pd.DataFrame:
    piv = stats.pivot(index="project", columns=label_col_name, values="delta_tumor_minus_normal")
    order = [o for o in order if o in piv.columns] if order else sorted(piv.columns)
    return piv[order]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--batch-summaries-glob",
        required=True,
        action="append",
        help="glob pattern matching per-batch summary CSVs; repeatable (e.g. one per processing run, such as a pilot cohort plus the full cohort)",
    )
    ap.add_argument("--region-types", default="AER,cts-ADR,ubi-ADR", help="comma-separated region-type categories to test (default: AER,cts-ADR,ubi-ADR)")
    ap.add_argument("--subgroup-order", default="", help="comma-separated preferred column order for the subgroup heatmap tables (default: discovery order)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    files = [f for pattern in args.batch_summaries_glob for f in glob.glob(pattern)]
    if not files:
        raise SystemExit("no files matched --batch-summaries-glob")
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    df = df.drop_duplicates(subset="sample_uuid")
    df.to_csv(outdir / "all_projects_summary.csv", index=False)
    print(f"total samples aggregated: {df.shape[0]}")

    region_cats = [r for r in args.region_types.split(",") if r]
    sg_cols_all = [c for c in df.columns if c.startswith("sg_") and c.endswith("_mean")]
    subgroups = [c[3:-5] for c in sg_cols_all]
    sg_order = [s for s in args.subgroup_order.split(",") if s] or subgroups

    # ---------- raw ----------
    stats = compute_stats(df, region_cats, lambda label: f"{label}_mean", "region_type")
    stats.to_csv(outdir / "all_projects_stats_regiontype.csv", index=False)
    sg_stats = compute_stats(df, subgroups, lambda label: f"sg_{label}_mean", "subgroup")
    sg_stats.to_csv(outdir / "all_projects_stats_subgroup.csv", index=False)

    heatmap_table(stats, "region_type", region_cats).to_csv(outdir / "heatmap_table_regiontype_delta.csv")
    heatmap_table(sg_stats, "subgroup", sg_order).to_csv(outdir / "heatmap_table_subgroup_delta.csv")

    # ---------- background (global-mean) normalized ----------
    for cat in region_cats:
        if f"{cat}_mean" in df.columns:
            df[f"{cat}_globalnorm"] = df[f"{cat}_mean"] - df["global_mean"]
    for sg in subgroups:
        if f"sg_{sg}_mean" in df.columns:
            df[f"sg_{sg}_globalnorm"] = df[f"sg_{sg}_mean"] - df["global_mean"]
    df.to_csv(outdir / "all_projects_summary_globalnorm.csv", index=False)

    global_stats = compute_stats(df, ["global"], lambda label: f"{label}_mean", "region_type")
    global_stats.to_csv(outdir / "all_projects_stats_global.csv", index=False)

    stats_gn = compute_stats(df, region_cats, lambda label: f"{label}_globalnorm", "region_type")
    stats_gn.to_csv(outdir / "all_projects_stats_regiontype_globalnorm.csv", index=False)
    sg_stats_gn = compute_stats(df, subgroups, lambda label: f"sg_{label}_globalnorm", "subgroup")
    sg_stats_gn.to_csv(outdir / "all_projects_stats_subgroup_globalnorm.csv", index=False)

    heatmap_table(stats_gn, "region_type", region_cats).to_csv(outdir / "heatmap_table_regiontype_delta_globalnorm.csv")
    heatmap_table(sg_stats_gn, "subgroup", sg_order).to_csv(outdir / "heatmap_table_subgroup_delta_globalnorm.csv")

    print("done")


if __name__ == "__main__":
    main()
