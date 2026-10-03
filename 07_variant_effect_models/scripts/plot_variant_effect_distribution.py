#!/usr/bin/env python
"""Distribution of predicted variant chromatin-accessibility effects.

Takes a variant x cell-type predicted delta-effect table (one `*_delta`
column per cell type) and characterizes the distribution of effects:

  1. pooled distribution of delta effects across all variants and cell
     types
  2. per-cell-type distribution (violin plots)
  3. per-variant cell-type specificity of delta effects (Gini index of
     |delta| across cell types -- 0 = effect spread evenly across cell
     types, 1 = concentrated in a single cell type)
  4. joint density of cell-type specificity (Gini) vs. effect magnitude
     (max |delta| across cell types)

Consolidated from seq2print/1-1j_enrichment_openTargets_study.ipynb
section 1.
"""

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def gini_per_row(x: np.ndarray) -> float:
    """Gini index of a non-negative vector (0 = perfectly even, 1 = fully concentrated)."""
    x = np.nan_to_num(np.asarray(x, dtype=float), nan=0.0)
    total = x.sum()
    if total == 0:
        return 0.0
    x = np.sort(x)
    n = len(x)
    return ((2 * np.arange(1, n + 1) - n - 1) @ x) / (n * total)


def load_long_delta(df_snp: pd.DataFrame, delta_cols: list[str]) -> pd.DataFrame:
    df_delta = df_snp[delta_cols].melt(var_name="cell_type", value_name="delta_effect")
    df_delta["cell_type"] = (
        df_delta["cell_type"].str.replace(r"^V2\.", "", regex=True).str.replace(r"_delta$", "", regex=True)
    )
    return df_delta


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--snp-delta", required=True, help="variant x celltype delta-effect table, tab-separated, with a `variant_id` column and one `*_delta` column per cell type")
    ap.add_argument("--celltype-colors", help="optional JSON mapping cell type name -> color, for the per-celltype violin plots (default: a single color for all)")
    ap.add_argument("--ncols", type=int, default=10, help="number of columns in the per-celltype violin grid (default: 10)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df_snp = pd.read_csv(args.snp_delta, sep="\t")
    delta_cols = [c for c in df_snp.columns if c.endswith("_delta")]
    if not delta_cols:
        raise SystemExit("no columns ending in '_delta' found in --snp-delta")
    print(f"{len(df_snp)} variants x {len(delta_cols)} cell types")

    df_delta = load_long_delta(df_snp, delta_cols)
    celltypes = df_delta["cell_type"].drop_duplicates().tolist()

    celltype_color_map = {}
    if args.celltype_colors:
        celltype_color_map = json.loads(Path(args.celltype_colors).read_text())
    celltype_colors = {ct: celltype_color_map.get(ct, "steelblue") for ct in celltypes}

    # ---------- 1-1 pooled distribution of delta effects ----------
    plt.figure(figsize=(9, 5))
    sns.histplot(data=df_delta, x="delta_effect", bins=40, kde=True, stat="density", color="steelblue")
    plt.axvline(0, color="black", linestyle="--", linewidth=1)
    plt.xlabel("Delta effect")
    plt.ylabel("Density")
    plt.title("Distribution of delta effects across variants and cell types")
    plt.tight_layout()
    plt.savefig(outdir / "distribution_delta_effect_pooled.png", dpi=300, bbox_inches="tight")
    plt.savefig(outdir / "distribution_delta_effect_pooled.pdf", bbox_inches="tight")
    plt.close()

    # ---------- 1-2 per-celltype distribution ----------
    ncols = args.ncols
    nrows = math.ceil(len(celltypes) / ncols)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(3 * ncols, 3.5 * nrows), sharex=True, sharey=True)
    axes = axes.flatten() if hasattr(axes, "flatten") else [axes]

    for ax, celltype in zip(axes, celltypes):
        values = df_delta.loc[df_delta["cell_type"] == celltype, "delta_effect"].dropna()
        sns.violinplot(y=values, ax=ax, color=celltype_colors[celltype], inner="box", cut=0)
        ax.axhline(0, color="black", linestyle="--", linewidth=0.8)
        ax.set_title(celltype.replace("_", " "))
        ax.set_xlabel("")
        ax.set_ylabel("Delta effect")
    for ax in axes[len(celltypes):]:
        ax.set_visible(False)

    fig.suptitle("Distribution of delta effects by cell type", fontsize=16, y=1.02)
    plt.tight_layout()
    fig.savefig(outdir / "distribution_delta_effect_per_celltype.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "distribution_delta_effect_per_celltype.pdf", bbox_inches="tight")
    plt.close(fig)

    # ---------- 1-3 celltype specificity (Gini index) ----------
    effect_matrix = df_snp[delta_cols].abs().to_numpy()
    df_snp["delta_gini"] = np.apply_along_axis(gini_per_row, axis=1, arr=effect_matrix)
    df_snp["n_nonzero_celltypes"] = (effect_matrix > 0).sum(axis=1)

    plt.figure(figsize=(8, 5))
    sns.histplot(data=df_snp, x="delta_gini", bins=100, kde=False, color="steelblue", edgecolor="white")
    plt.axvline(
        df_snp["delta_gini"].mean(), color="darkred", linestyle="--", linewidth=1.5,
        label=f"Mean = {df_snp['delta_gini'].mean():.3f}",
    )
    plt.xlabel("Gini index of absolute delta effects")
    plt.ylabel("Number of variants")
    plt.title("Cell-type specificity of delta effects across variants")
    plt.legend(frameon=False)
    plt.grid(False)
    sns.despine()
    plt.tight_layout()
    plt.savefig(outdir / "distribution_delta_effect_gini.png", dpi=300, bbox_inches="tight")
    plt.savefig(outdir / "distribution_delta_effect_gini.pdf", bbox_inches="tight")
    plt.close()

    # ---------- 1-4 Gini vs. max |delta| ----------
    df_snp["max_abs_delta"] = df_snp[delta_cols].abs().max(axis=1)
    plot_df = df_snp[["variant_id", "delta_gini", "max_abs_delta"]].replace([np.inf, -np.inf], np.nan).dropna()

    sns.set_style("white")
    fig, ax = plt.subplots(figsize=(7, 6))
    hb = ax.hexbin(plot_df["delta_gini"], plot_df["max_abs_delta"], gridsize=35, mincnt=1, bins="log", cmap="viridis")
    cbar = fig.colorbar(hb, ax=ax)
    cbar.set_label("Number of variants per bin")
    ax.set_xlabel("Gini index of absolute delta effects")
    ax.set_ylabel("Maximum |delta effect|")
    ax.set_title("Density of variants by cell-type specificity and maximum effect")
    ax.grid(False)
    sns.despine()
    plt.tight_layout()
    plt.savefig(outdir / "gini_vs_max_delta.png", dpi=300, bbox_inches="tight")
    plt.savefig(outdir / "gini_vs_max_delta.pdf", bbox_inches="tight")
    plt.close(fig)

    df_snp[["variant_id", "delta_gini", "n_nonzero_celltypes", "max_abs_delta"]].to_csv(
        outdir / "variant_effect_specificity.tsv", sep="\t", index=False
    )

    print("done")


if __name__ == "__main__":
    main()
