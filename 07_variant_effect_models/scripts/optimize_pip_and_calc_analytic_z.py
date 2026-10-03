#!/usr/bin/env python
"""Optimize the GWAS fine-mapping PIP threshold per study and compute an
analytic z-test of seq2PRINT variant effects per (study, cell type).

For each GWAS study and each candidate PIP threshold, computes the mean
|predicted delta| across that study's credible-set SNPs passing the
threshold, per cell type (`compute_stats_multi_pip`). Each (study,
threshold, celltype) mean is then analytic-z-tested against a background
built from every *other* study's SNPs at that (threshold, celltype) --
i.e. a leave-one-study-out null (`analytic_z_test`). Because the "right"
PIP threshold trades off SNP count against fine-mapping confidence and
isn't known a priori, the threshold that gives each study its clearest
enrichment (highest median z across cell types) is picked as that study's
own optimized threshold (`pick_best_pip_threshold`), and the final table
is filtered down to each study's one optimized-threshold row per cell
type (`filter_to_best_threshold`).

Two additional checks from the notebook are available but off by default
(they roughly double runtime):
  --declump-sensitivity-check: repeats the whole threshold-sweep +
    z-test on a LD-decorrelated SNP set (one lead SNP per ~500kb locus
    per study), to check whether enrichment is an LD-driven artifact of
    treating correlated credible-set SNPs as independent observations.
  --sanity-check-distributions: renders a multi-page PDF of each study's
    actual (signed) per-celltype delta distribution at its optimized
    threshold, to flag means driven by a few outlier SNPs.

Requires per-study SNP-level delta tables already built (one row per
credible-set SNP, columns PIP, chrom, pos, and one `*_delta` column per
cell type) -- joining GWAS fine-mapping output with seq2PRINT predicted
variant effects is a separate, heavier data-prep step and out of scope
here.

Consolidated from seq2print/1-1f-1_permutate_optimize_seq2print_variant.ipynb
section 1 (1-1 through 1-5).
"""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import norm
from tqdm import tqdm

LOCUS_BIN_SIZE = 500_000  # bp; crude proxy for one independent LD block


def declump_lead_snps(snp_df: pd.DataFrame, bin_size: int = LOCUS_BIN_SIZE) -> pd.DataFrame:
    """Collapse to one lead (max-PIP) SNP per ~bin_size-bp window per chromosome."""
    df = snp_df.copy()
    df["_locus"] = df["chrom"].astype(str) + ":" + (df["pos"] // bin_size).astype(str)
    lead = df.sort_values("PIP", ascending=False).drop_duplicates(subset="_locus", keep="first")
    return lead.drop(columns="_locus")


def compute_stats_multi_pip(df_efo, delta_dir: Path, pip_thresholds, min_snps: int, declump: bool = False):
    """Mean |delta| per (study, celltype) across a sweep of PIP thresholds.

    Returns (combined_stats, per_study_deltas_by_thr): combined_stats has
    one row per (study, celltype, pip_threshold); per_study_deltas_by_thr
    maps pip_threshold -> celltype -> {study_id: |delta| array}, the pool
    used to build each study's leave-one-study-out background.
    """
    all_stats_list = []
    per_study_deltas_by_thr = {}

    for pip_thr in pip_thresholds:
        per_study_deltas = {}
        per_thr_stats = []

        for _, row in df_efo.iterrows():
            study_id = row["study_id"]
            efo_id = row["efo_id"]
            phenotype = row["phenotype"]
            safe_pheno = str(phenotype).replace(" ", "_").replace("/", "_")

            snp_path = delta_dir / f"snp_with_delta.{study_id}.tsv"
            if not snp_path.exists():
                continue

            snp_df = pd.read_csv(snp_path, sep="\t")
            if declump:
                snp_df = declump_lead_snps(snp_df)

            snp_df_thr = snp_df[snp_df["PIP"] >= pip_thr].copy()
            n_snps = snp_df_thr.shape[0]
            if n_snps <= min_snps:
                continue

            delta_cols = [c for c in snp_df_thr.columns if c.endswith("_delta")]
            if not delta_cols:
                continue

            delta_abs = np.abs(snp_df_thr[delta_cols])
            for col in delta_cols:
                ct = col.replace("_delta", "").replace("V2.", "")
                per_study_deltas.setdefault(ct, {})[study_id] = delta_abs[col].values

            col_mean = delta_abs.mean(axis=0)
            stats = pd.DataFrame({
                "celltype": [c.replace("_delta", "").replace("V2.", "") for c in delta_cols],
                "mean_abs_delta": col_mean.values,
            })
            stats["study_id"] = study_id
            stats["efo_id"] = efo_id
            stats["pheno"] = safe_pheno
            stats["n_snps"] = n_snps
            stats["pip_threshold"] = pip_thr
            per_thr_stats.append(stats)

        per_thr_stats_df = (
            pd.concat(per_thr_stats, ignore_index=True)
            if per_thr_stats
            else pd.DataFrame(columns=["celltype", "mean_abs_delta", "study_id", "efo_id", "pheno", "n_snps", "pip_threshold"])
        )
        per_study_deltas_by_thr[pip_thr] = per_study_deltas
        all_stats_list.append(per_thr_stats_df)

    combined_stats = pd.concat(all_stats_list, ignore_index=True)
    return combined_stats, per_study_deltas_by_thr


def _get_background(per_study_deltas_by_thr, pip_thr, celltype, exclude_study_id):
    per_study = per_study_deltas_by_thr.get(pip_thr, {}).get(celltype, {})
    arrays = [arr for sid, arr in per_study.items() if sid != exclude_study_id]
    return np.concatenate(arrays) if arrays else np.array([])


def analytic_z_test(combined_stats: pd.DataFrame, per_study_deltas_by_thr: dict, min_snps: int) -> pd.DataFrame:
    """Analytic z-test of each (study, celltype, threshold)'s mean |delta| against a
    leave-one-study-out background pooled from every other study."""
    combined_stats = combined_stats.copy()
    null_mean_list, null_sd_list, z_list, p_list = [], [], [], []

    for _, row in tqdm(combined_stats.iterrows(), total=len(combined_stats), desc="Analytic z-test"):
        pip_thr, celltype, study_id = row["pip_threshold"], row["celltype"], row["study_id"]
        n_snps, obs_mean = int(row["n_snps"]), float(row["mean_abs_delta"])

        if n_snps < min_snps:
            null_mean_list.append(np.nan); null_sd_list.append(np.nan); z_list.append(np.nan); p_list.append(np.nan)
            continue

        bg = _get_background(per_study_deltas_by_thr, pip_thr, celltype, exclude_study_id=study_id)
        if len(bg) == 0:
            null_mean_list.append(np.nan); null_sd_list.append(np.nan); z_list.append(np.nan); p_list.append(np.nan)
            continue

        mu0 = bg.mean()
        sd_bg = bg.std(ddof=1)
        se = sd_bg / np.sqrt(n_snps)
        if se > 0:
            z = (obs_mean - mu0) / se
            p = 2 * (1 - norm.cdf(abs(z)))
        else:
            z, p = np.nan, np.nan

        null_mean_list.append(mu0); null_sd_list.append(sd_bg); z_list.append(z); p_list.append(p)

    combined_stats["null_mean_analytic"] = null_mean_list
    combined_stats["null_sd_analytic"] = null_sd_list
    combined_stats["z_analytic"] = z_list
    combined_stats["p_analytic_two_sided"] = p_list
    return combined_stats


def compute_study_threshold_scores(combined_stats, study_id, min_snps, use_common_celltypes=True):
    df = combined_stats[
        (combined_stats["study_id"] == study_id)
        & (combined_stats["n_snps"] >= min_snps)
        & (~combined_stats["z_analytic"].isna())
    ].copy()
    if df.empty:
        return None

    z_by_thr = {thr: g.set_index("celltype")["z_analytic"] for thr, g in df.groupby("pip_threshold")}
    thresholds = sorted(z_by_thr.keys())

    if use_common_celltypes:
        common = set.intersection(*[set(s.index) for s in z_by_thr.values()])
        if len(common) >= 3:
            for thr in thresholds:
                z_by_thr[thr] = z_by_thr[thr].loc[list(common)]

    rows = []
    for thr in thresholds:
        zvals = z_by_thr[thr].values
        rows.append({
            "study_id": study_id, "pip_threshold": thr,
            "score_median_z": np.median(zvals), "n_celltypes_used": len(zvals),
        })
    return pd.DataFrame(rows).sort_values("pip_threshold")


def pick_best_pip_threshold(combined_stats: pd.DataFrame, min_snps: int):
    """Per study, the PIP threshold with the highest median z_analytic across
    cell types shared by every threshold ("clearest enrichment")."""
    study_ids = combined_stats["study_id"].unique()
    score_curves, best_rows = [], []

    for sid in tqdm(study_ids, desc="Picking best PIP threshold per study"):
        sdf = compute_study_threshold_scores(combined_stats, sid, min_snps, use_common_celltypes=True)
        if sdf is None or sdf.empty:
            continue
        score_curves.append(sdf)
        best = sdf.loc[sdf["score_median_z"].idxmax()]
        best_rows.append({
            "study_id": sid,
            "best_pip_threshold": float(best["pip_threshold"]),
            "best_score_median_z": float(best["score_median_z"]),
            "n_celltypes_used_at_best": int(best["n_celltypes_used"]),
        })

    score_curves_df = pd.concat(score_curves, ignore_index=True) if score_curves else pd.DataFrame()
    best_per_study_df = pd.DataFrame(best_rows)
    return score_curves_df, best_per_study_df


def filter_to_best_threshold(combined_stats: pd.DataFrame, best_per_study_df: pd.DataFrame) -> pd.DataFrame:
    best_map = best_per_study_df[["study_id", "best_pip_threshold"]].copy()
    cs_best = combined_stats.merge(best_map, on="study_id", how="inner")
    cs_best = cs_best[np.isclose(cs_best["pip_threshold"], cs_best["best_pip_threshold"])].copy()
    return cs_best


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--efo-studies", required=True, help="tab-separated table of GWAS studies to process: columns study_id, efo_id, phenotype")
    ap.add_argument("--delta-perstudy-dir", required=True, help="directory of per-study SNP-level delta tables, named snp_with_delta.{study_id}.tsv")
    ap.add_argument("--pip-thresholds", default="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0", help="comma-separated PIP thresholds to sweep (default: 0.1..1.0 by 0.1)")
    ap.add_argument("--min-snps", type=int, default=50, help="minimum SNPs required to test a (study, threshold) (default: 50)")
    ap.add_argument("--declump-sensitivity-check", action="store_true", help="also repeat the analysis on one lead SNP per ~500kb locus per study, and plot it against the main z_analytic (default: off, roughly doubles runtime)")
    ap.add_argument("--sanity-check-distributions", action="store_true", help="also render a multi-page PDF of each study's per-celltype delta distribution at its optimized threshold (default: off)")
    ap.add_argument("--celltype-colors", help="optional JSON mapping cell type name -> color, for --sanity-check-distributions (default: gray)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    delta_dir = Path(args.delta_perstudy_dir)
    pip_thresholds = [float(x) for x in args.pip_thresholds.split(",")]
    df_efo = pd.read_csv(args.efo_studies, sep="\t")

    # ---------- 1-1 / 1-2 ----------
    print("1-1: computing mean |delta| across PIP thresholds")
    combined_stats, per_study_deltas_by_thr = compute_stats_multi_pip(df_efo, delta_dir, pip_thresholds, args.min_snps)
    combined_stats.to_csv(outdir / "all_SNP_delta_stats_multiPIP.tsv", sep="\t", index=False)

    print("1-2: analytic z-test vs. leave-one-study-out background")
    combined_stats = analytic_z_test(combined_stats, per_study_deltas_by_thr, args.min_snps)
    combined_stats.to_csv(outdir / "all_SNP_delta_stats_multiPIP.with_analytic_z.tsv", sep="\t", index=False)

    # ---------- 1-2b ----------
    if args.declump_sensitivity_check:
        print("1-2b: LD-declumped sensitivity check")
        combined_stats_locus, per_study_deltas_by_thr_locus = compute_stats_multi_pip(
            df_efo, delta_dir, pip_thresholds, args.min_snps, declump=True
        )
        combined_stats_locus = analytic_z_test(combined_stats_locus, per_study_deltas_by_thr_locus, args.min_snps)
        combined_stats_locus.to_csv(outdir / "all_SNP_delta_stats_multiPIP.locus_declumped.with_analytic_z.tsv", sep="\t", index=False)

        merged_locus = combined_stats.merge(
            combined_stats_locus[["study_id", "celltype", "pip_threshold", "z_analytic"]].rename(
                columns={"z_analytic": "z_analytic_locus"}
            ),
            on=["study_id", "celltype", "pip_threshold"], how="inner",
        )
        r = merged_locus[["z_analytic", "z_analytic_locus"]].corr().iloc[0, 1]

        plt.figure(figsize=(5, 5))
        plt.scatter(merged_locus["z_analytic"], merged_locus["z_analytic_locus"], s=8, alpha=0.4)
        lims = [merged_locus[["z_analytic", "z_analytic_locus"]].min().min(), merged_locus[["z_analytic", "z_analytic_locus"]].max().max()]
        plt.plot(lims, lims, color="grey", ls="--")
        plt.xlabel("z_analytic (all SNPs, 1-2)")
        plt.ylabel("z_analytic (one SNP per locus, 1-2b)")
        plt.title(f"Robustness to LD (r = {r:.2f})")
        plt.tight_layout()
        plt.savefig(figdir / "z_analytic_LD_robustness.png", dpi=300)
        plt.savefig(figdir / "z_analytic_LD_robustness.pdf")
        plt.close()

    # ---------- 1-3 ----------
    print("1-3: picking the best PIP threshold per study")
    score_curves_df, best_per_study_df = pick_best_pip_threshold(combined_stats, args.min_snps)
    score_curves_df.to_csv(outdir / "study_threshold_score_curves.tsv", sep="\t", index=False)
    best_per_study_df.to_csv(outdir / "best_pip_threshold_per_study.tsv", sep="\t", index=False)

    plt.figure()
    sns.histplot(data=best_per_study_df, x="best_pip_threshold", bins=30, stat="count", color="gray", edgecolor="black")
    plt.xlim(0.1, 1.0)
    plt.xticks(np.arange(0.1, 1.01, 0.1))
    plt.xlabel("Optimized PIP threshold")
    plt.ylabel("Number of study")
    plt.title("Distribution of optimized PIP thresholds per study")
    plt.tight_layout()
    plt.savefig(figdir / "optimized_pip_thresholds_hist.png", dpi=300)
    plt.savefig(figdir / "optimized_pip_thresholds_hist.pdf")
    plt.close()

    # ---------- 1-4 ----------
    print("1-4: filtering to each study's optimized threshold")
    cs_best = filter_to_best_threshold(combined_stats, best_per_study_df)
    cs_best.to_csv(outdir / "delta_score.filtered_to_best_threshold.tsv", sep="\t", index=False)
    print(f"{cs_best.shape[0]} rows, {cs_best['study_id'].nunique()} studies")

    # ---------- 1-5 ----------
    if args.sanity_check_distributions:
        print("1-5: sanity-check distributions")
        celltype_color_map = json.loads(Path(args.celltype_colors).read_text()) if args.celltype_colors else {}
        pheno_map = cs_best.groupby("study_id")["pheno"].first().to_dict()
        best_map = best_per_study_df[["study_id", "best_pip_threshold"]].copy()

        from matplotlib.backends.backend_pdf import PdfPages
        from matplotlib.lines import Line2D

        out_pdf = figdir / "seq2PRINT_delta_distribution_per_study.bestPIP.pdf"
        with PdfPages(out_pdf) as pdf:
            legend_cts = sorted(celltype_color_map.keys())
            if legend_cts:
                handles = [Line2D([0], [0], color=celltype_color_map[ct], lw=2) for ct in legend_cts]
                plt.figure(figsize=(6, 0.22 * len(legend_cts) + 1))
                plt.legend(handles, legend_cts, loc="center", ncol=2, fontsize=7, frameon=False)
                plt.axis("off")
                plt.title("Celltype color legend")
                pdf.savefig()
                plt.close()

            for _, row in tqdm(best_map.iterrows(), total=len(best_map), desc="Writing per-study delta distributions"):
                study_id, pip_thr = row["study_id"], row["best_pip_threshold"]
                snp_path = delta_dir / f"snp_with_delta.{study_id}.tsv"
                if not snp_path.exists():
                    continue
                snp_df = pd.read_csv(snp_path, sep="\t")
                snp_df_thr = snp_df[snp_df["PIP"] >= pip_thr].copy()
                if snp_df_thr.empty:
                    continue
                delta_cols = [c for c in snp_df_thr.columns if c.endswith("_delta")]
                if not delta_cols:
                    continue

                df_long = snp_df_thr[delta_cols].melt(var_name="celltype", value_name="delta")
                df_long["celltype"] = df_long["celltype"].str.replace("_delta", "", regex=False).str.replace("V2.", "", regex=False)
                pheno = pheno_map.get(study_id, "")

                plt.figure(figsize=(10, 8))
                for ct, sub in df_long.groupby("celltype"):
                    color = celltype_color_map.get(ct, "gray")
                    sns.kdeplot(sub["delta"], color=color, linewidth=1.2, alpha=0.85)
                plt.axvline(0, color="black", linestyle="--", linewidth=1)
                plt.xlabel("seq2PRINT delta (signed)")
                plt.ylabel("Density")
                plt.title(f"{study_id}\n{pheno}\nn_snps={snp_df_thr.shape[0]}, PIP>={pip_thr}")
                plt.tight_layout()
                pdf.savefig()
                plt.close()
        print(f"saved {out_pdf}")

    print("done")


if __name__ == "__main__":
    main()
