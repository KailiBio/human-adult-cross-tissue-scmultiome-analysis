#!/usr/bin/env python
"""Volcano plots of TF-expression-vs-ChromVAR correlation.

Takes the per-(TF, celltype) correlation table produced by
compute_tf_chromvar_correlation.py (the per-celltype-FDR version, i.e.
`tf_chromvar_correlation.per_celltype.ctFDR.txt`) and makes two kinds of
volcano plot:

  - one per TF, across cell types (r vs. -log10(FDR), one point per
    celltype)
  - one per cell type, across TFs (r vs. -log10(FDR), one point per TF)

Consolidated from seq2print/3-2_ChromVAR_TFexpression_correlation.ipynb
section 1 (1-3-1 volcano_per_tf_celltype, 1-3-2 loop, 1-3-3
volcano_per_celltype_all_tfs + loop).
"""

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common.naming import sanitize_label


def volcano_per_tf_celltype(
    df_r2,
    tf_name,
    r_col="r",
    r2_col="r2",
    p_col="fdr",
    top_n=15,
    palette=None,
    outdir=".",
):
    """Volcano for ONE TF, across cell types."""
    df = df_r2.dropna(subset=[r_col, r2_col, p_col]).copy()
    if df.empty:
        print(f"{tf_name}: no rows to plot.")
        return

    df["minus_log10_p"] = -np.log10(df[p_col].where(df[p_col] > 0))

    top_by_p = df.nlargest(top_n, "minus_log10_p")
    top_by_r2 = df.nlargest(top_n, r2_col)
    top_union = pd.concat([top_by_p, top_by_r2]).drop_duplicates(subset=["celltype"])
    df["is_top"] = df["celltype"].isin(top_union["celltype"])

    plt.figure(figsize=(10, 6))
    ax = sns.scatterplot(
        data=df, x=r_col, y="minus_log10_p", hue="celltype", palette=palette,
        s=45, alpha=0.9, edgecolor="black", linewidth=0.3, legend=False,
    )

    a = round(df[r_col].abs().max(), 1) + 0.1
    ax.set_xlim((-a, a) if a < 0.8 else (-0.8, 0.8))

    plt.axvline(0, color="grey", linestyle="--", linewidth=1)
    for _, row in df[df["is_top"]].iterrows():
        plt.text(row[r_col], row["minus_log10_p"], row["celltype"], fontsize=8, ha="center", va="bottom")

    plt.xlabel("Spearman r (expr vs chromVAR, per cell type)")
    plt.ylabel("-log10(FDR)")
    plt.title(tf_name)
    plt.grid(False)
    plt.tight_layout()

    plt.savefig(f"{outdir}/TF_volcano_by_celltype.{tf_name}.volcano.png", dpi=300)
    plt.savefig(f"{outdir}/TF_volcano_by_celltype.{tf_name}.volcano.pdf", dpi=300)
    plt.close()


def volcano_per_celltype_all_tfs(
    df_all,
    celltype,
    r_col="r",
    r2_col="r2",
    p_col="fdr",
    top_n=15,
    outdir=".",
):
    """Volcano for ONE cell type, across all TFs."""
    df = df_all[df_all["celltype"] == celltype].copy()
    if df.empty:
        print(f"No rows for cell type: {celltype}")
        return

    df = df.dropna(subset=[r_col, r2_col, p_col]).copy()
    if df.empty:
        print(f"No valid r / p rows for cell type: {celltype}")
        return

    df["minus_log10_p"] = -np.log10(df[p_col].where(df[p_col] > 0))

    top_by_p = df.nlargest(top_n, "minus_log10_p")
    top_by_r2 = df.nlargest(top_n, r2_col)
    top_union = pd.concat([top_by_p, top_by_r2]).drop_duplicates(subset=["TF"])
    df["is_top"] = df["TF"].isin(top_union["TF"])

    plt.figure(figsize=(10, 7))
    ax = sns.scatterplot(
        data=df, x=r_col, y="minus_log10_p", s=35, alpha=0.8,
        edgecolor="black", linewidth=0.3, color="steelblue", legend=False,
    )
    ax.set_xlim(-1, 1)
    ax.set_xticks([-1, -0.8, -0.5, -0.3, 0, 0.3, 0.5, 0.8, 1])

    plt.axvline(0, color="grey", linestyle="--", linewidth=1)
    for _, row in df[df["is_top"]].iterrows():
        plt.text(row[r_col], row["minus_log10_p"], row["TF"], fontsize=11, ha="center", va="bottom")

    plt.xlabel("Spearman r (expr vs chromVAR, per TF)")
    plt.ylabel("-log10(FDR)")
    plt.title(f"{celltype} (n={df['n_cells'].iloc[0]})")
    plt.grid(False)
    plt.tight_layout()

    safe_ct = sanitize_label(celltype)
    plt.savefig(f"{outdir}/volcano_allTFs.celltype_{safe_ct}.png", dpi=300)
    plt.savefig(f"{outdir}/volcano_allTFs.celltype_{safe_ct}.pdf", dpi=300)
    plt.close()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--corr-table",
        required=True,
        help="per-celltype-FDR correlation table, e.g. tf_chromvar_correlation.per_celltype.ctFDR.txt "
        "from compute_tf_chromvar_correlation.py (columns: celltype, TF, r, r2, p, n_cells, fdr, ...)",
    )
    ap.add_argument(
        "--celltype-colors",
        help="optional tab-separated table (celltype, color) to color points by celltype in the per-TF volcano plots",
    )
    ap.add_argument("--min-cells-per-tf-volcano", type=int, default=100, help="minimum n_cells for a (TF, celltype) point in the per-TF volcano plots (default: 100)")
    ap.add_argument("--top-n", type=int, default=15, help="number of top points (by -log10(FDR) and by R^2) to label (default: 15)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    tf_volcano_dir = outdir / "TF_volcano"
    celltype_volcano_dir = outdir / "celltype_volcano"
    tf_volcano_dir.mkdir(parents=True, exist_ok=True)
    celltype_volcano_dir.mkdir(parents=True, exist_ok=True)

    df_all = pd.read_csv(args.corr_table, sep="\t")

    palette = None
    if args.celltype_colors:
        colors = pd.read_csv(args.celltype_colors, sep="\t", header=None, names=["celltype", "color"])
        palette = dict(zip(colors["celltype"], colors["color"]))

    tfs = df_all["TF"].unique()
    print(f"plotting per-TF volcanoes for {len(tfs)} TFs")
    for tf in tfs:
        df_tf = df_all[df_all["TF"] == tf]
        df_big = df_tf[df_tf["n_cells"] >= args.min_cells_per_tf_volcano]
        volcano_per_tf_celltype(df_big, tf_name=tf, top_n=args.top_n, palette=palette, outdir=str(tf_volcano_dir))

    celltypes = sorted(df_all["celltype"].unique())
    print(f"plotting per-celltype volcanoes for {len(celltypes)} cell types")
    for ct in celltypes:
        volcano_per_celltype_all_tfs(df_all, celltype=ct, top_n=args.top_n, outdir=str(celltype_volcano_dir))

    print("done")


if __name__ == "__main__":
    main()
