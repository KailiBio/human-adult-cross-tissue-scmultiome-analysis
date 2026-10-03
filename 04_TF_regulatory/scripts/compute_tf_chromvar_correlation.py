#!/usr/bin/env python
"""Correlate TF expression (RNA) with its own ChromVAR motif deviation score.

For each shared TF/gene and each cell type, computes the Spearman
correlation (r, R^2, p) between that TF's RNA expression and its ChromVAR
deviation score across cells of that type (with per-celltype downsampling,
to keep large cell types from dominating runtime). Results are then FDR
(Benjamini-Hochberg) corrected two ways: within each celltype (across TFs)
and globally (across every TF x celltype instance).

Expects `--rna-h5ad`/`--chromvar-h5ad` to already share the same cell
barcodes (obs_names) and the same TF/gene set (var_names) -- reconciling
RNA cell ids with ATAC/ChromVAR barcodes is a separate, cohort-specific
data-prep step done upstream of this script.

Consolidated from seq2print/3-2_ChromVAR_TFexpression_correlation.ipynb
section 1 (1-3-1, 1-3-2 compute half; FDR correction cells).
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


def to_dense(X):
    return X.toarray() if hasattr(X, "toarray") else X


def compute_tf_r2_per_celltype(
    tf_name,
    rna_ad,
    chromvar_ad,
    celltype_key="CellAnnotation_L1",
    min_cells=20,
    max_cells_per_celltype=None,
    random_state=0,
    outdir=None,
):
    """Spearman r/R^2 between TF expression and ChromVAR, within each cell type.

    Returns a DataFrame with columns: celltype, r, r2, p, n_cells, TF.
    """
    if tf_name not in rna_ad.var_names or tf_name not in chromvar_ad.var_names:
        raise ValueError(f"{tf_name} not found in both RNA and chromVAR var_names")

    x_expr = to_dense(rna_ad[:, tf_name].X).ravel()
    y_chromvar = to_dense(chromvar_ad[:, tf_name].X).ravel()
    celltypes = rna_ad.obs[celltype_key].astype(str).values

    df_cells = pd.DataFrame({"expression": x_expr, "chromVAR": y_chromvar, "celltype": celltypes})

    empty_cols = ["celltype", "r", "r2", "p", "n_cells", "TF"]
    if df_cells.empty:
        print(f"{tf_name}: no cells found to analyze.")
        empty = pd.DataFrame(columns=empty_cols)
        if outdir:
            empty.to_csv(f"{outdir}/TF_R_per_celltype.{tf_name}.txt", sep="\t", index=False)
        return empty

    rng = np.random.default_rng(random_state)
    results = []

    for ct, df_ct in df_cells.groupby("celltype"):
        n_total = df_ct.shape[0]
        if n_total < min_cells:
            continue

        if max_cells_per_celltype is not None and n_total > max_cells_per_celltype:
            idx = rng.choice(df_ct.index.values, size=max_cells_per_celltype, replace=False)
            df_ct = df_ct.loc[idx]
        n = df_ct.shape[0]

        x = df_ct["expression"].values
        y = df_ct["chromVAR"].values

        if np.std(x) == 0 or np.std(y) == 0:
            r, p = 0.0, np.nan
        else:
            r, p = spearmanr(x, y)

        results.append({"celltype": ct, "r": r, "r2": r**2, "p": p, "n_cells": n})

    df_r2 = pd.DataFrame(results)
    if df_r2.empty:
        print(f"{tf_name}: no cell types passed filtering / variance checks.")
        if outdir:
            df_r2.to_csv(f"{outdir}/TF_R_per_celltype.{tf_name}.txt", sep="\t", index=False)
        return df_r2

    df_r2["TF"] = tf_name
    df_r2 = df_r2.sort_values("r2", ascending=False)

    if outdir:
        df_r2.to_csv(f"{outdir}/TF_R_per_celltype.{tf_name}.txt", sep="\t", index=False)

    return df_r2


def add_fdr_per_celltype(df_all: pd.DataFrame, alpha: float) -> pd.DataFrame:
    """BH-FDR correct p-values within each celltype (across TFs)."""
    eps = 1e-300
    df_list = []
    for ct, df_ct in df_all.groupby("celltype"):
        df_ct = df_ct.copy()
        p = df_ct["p"].values
        mask = np.isfinite(p) & (p >= 0) & (p <= 1)

        df_ct["fdr"] = np.nan
        df_ct["minus_log10_FDR"] = np.nan
        df_ct["signed_minus_log10_FDR"] = np.nan
        df_ct["signif_fdr_ct"] = False

        if mask.sum() > 0:
            rej, p_fdr_ct, _, _ = multipletests(p[mask], alpha=alpha, method="fdr_bh")
            df_ct.loc[mask, "fdr"] = p_fdr_ct
            df_ct.loc[mask, "signif_fdr_ct"] = rej

            minus_log10 = -np.log10(np.clip(p_fdr_ct, eps, 1.0))
            df_ct.loc[mask, "minus_log10_FDR"] = minus_log10
            df_ct.loc[mask, "signed_minus_log10_FDR"] = np.sign(df_ct.loc[mask, "r"].values) * minus_log10

        df_list.append(df_ct)

    return pd.concat(df_list, ignore_index=True)


def add_fdr_global(df_all: pd.DataFrame, alpha: float) -> pd.DataFrame:
    """BH-FDR correct p-values globally, across every TF x celltype instance."""
    eps = 1e-300
    df_all = df_all.copy()
    pvals = df_all["p"].values
    mask = np.isfinite(pvals) & (pvals >= 0) & (pvals <= 1)

    rej, pvals_fdr_valid, _, _ = multipletests(pvals[mask], alpha=alpha, method="fdr_bh")
    pvals_fdr = np.full_like(pvals, np.nan, dtype=float)
    pvals_fdr[mask] = pvals_fdr_valid

    df_all["fdr"] = pvals_fdr
    df_all["minus_log10_FDR"] = -np.log10(np.clip(df_all["fdr"], eps, 1.0))
    df_all["signed_minus_log10_FDR"] = np.sign(df_all["r"]) * df_all["minus_log10_FDR"]

    return df_all


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rna-h5ad", required=True, help="RNA AnnData, cells x genes, aligned to --chromvar-h5ad (same obs_names)")
    ap.add_argument("--chromvar-h5ad", required=True, help="ChromVAR AnnData, cells x TFs, aligned to --rna-h5ad (same obs_names)")
    ap.add_argument("--celltype-key", default="CellAnnotation_L1", help="obs column with cell type labels (default: CellAnnotation_L1)")
    ap.add_argument("--min-cells", type=int, default=20, help="minimum cells required to correlate a TF within a celltype (default: 20)")
    ap.add_argument("--max-cells-per-celltype", type=int, default=5000, help="downsample each celltype to at most this many cells before correlating (default: 5000)")
    ap.add_argument("--random-state", type=int, default=1, help="seed for downsampling (default: 1)")
    ap.add_argument("--fdr-alpha", type=float, default=0.01, help="FDR (BH) alpha (default: 0.01)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    per_tf_dir = outdir / "per_tf"
    per_tf_dir.mkdir(parents=True, exist_ok=True)

    rna_ad = sc.read_h5ad(args.rna_h5ad)
    chromvar_ad = sc.read_h5ad(args.chromvar_h5ad)
    shared_tfs = chromvar_ad.var_names.tolist()

    print(f"correlating {len(shared_tfs)} TFs against ChromVAR across '{args.celltype_key}' groups")
    all_results = []
    for tf in shared_tfs:
        df_tf = compute_tf_r2_per_celltype(
            tf_name=tf,
            rna_ad=rna_ad,
            chromvar_ad=chromvar_ad,
            celltype_key=args.celltype_key,
            min_cells=args.min_cells,
            max_cells_per_celltype=args.max_cells_per_celltype,
            random_state=args.random_state,
            outdir=str(per_tf_dir),
        )
        all_results.append(df_tf)

    df_all = pd.concat(all_results, ignore_index=True)
    df_all.to_csv(outdir / "tf_chromvar_correlation.per_celltype.txt", sep="\t", index=False)

    df_ctfdr = add_fdr_per_celltype(df_all, args.fdr_alpha)
    df_ctfdr.to_csv(outdir / "tf_chromvar_correlation.per_celltype.ctFDR.txt", sep="\t", index=False)

    df_global = add_fdr_global(df_all, args.fdr_alpha)
    df_global.to_csv(outdir / "tf_chromvar_correlation.per_celltype.globalFDR.txt", sep="\t", index=False)

    print(f"done: {len(df_all)} TF x celltype instances")


if __name__ == "__main__":
    main()
