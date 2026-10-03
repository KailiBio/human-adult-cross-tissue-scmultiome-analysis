#!/usr/bin/env python
"""Celltype x tissue composition figures: cell counts, composition heatmaps,
and ranked tissue-specificity metrics.

From a per-cell table of (tissue, celltype) labels, produces:

  1. a horizontal bar chart of cells per celltype (log10, largest on top)
  2. two composition heatmaps: celltype breakdown of each tissue (rows sum
     to 100% per tissue) and tissue breakdown of each celltype (columns sum
     to 100% per celltype)
  3. four ranked scatter plots of each celltype's distribution across
     tissues, computed from the tissue-normalized (row-wise) composition
     table: its single highest per-tissue percentage ("max %"), its third
     highest ("3rd max %"), its median percentage across tissues, and the
     Gini index of that percentage across tissues (0 = evenly spread,
     1 = concentrated in one tissue). Median and Gini additionally get a
     labeled ("every point annotated with its celltype") version.
  4. celltype_composition_metrics.tsv, the four per-celltype metrics above
     in one table -- not written out by the source notebook (only plotted),
     added here since the values are already computed and a table is more
     reusable downstream than a figure.

Note on cell 158's original heatmap in the source notebook: it plotted
proportions_col.T (celltype-rows x tissue-columns) under axis labels and an
inline comment that both describe the UN-transposed orientation (tissue-rows
x celltype-columns) -- a copy-paste leftover from the preceding heatmap cell.
This script uses the un-transposed orientation, matching the labels/comment.

Consolidated from 0-9-2_concatenate_and_make_UMAP.ipynb section 3-2-1.
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common.colors import default_color_dict

DEFAULT_HEATMAP_COLORS = [
    "#000000", "#22350F", "#3B680C", "#529111", "#6ABB15",
    "#74DD0F", "#7FFF00", "#ADF121", "#D1E131", "#FFFFFF",
]


def gini_index(values) -> float:
    """Gini index of a non-negative vector (0 = perfectly even, 1 = fully concentrated)."""
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0 or np.sum(x) == 0:
        return np.nan
    if np.any(x < 0):
        raise ValueError("Gini index requires non-negative values.")
    x = np.sort(x)
    n = len(x)
    return (2 * np.sum(np.arange(1, n + 1) * x) / (n * np.sum(x))) - (n + 1) / n


def compute_counts_and_proportions(df: pd.DataFrame, tissue_col: str, celltype_col: str):
    """celltype_counts: overall cells per celltype, ascending (so a horizontal
    bar chart puts the largest celltype on top). counts_ct: tissue x celltype
    count table, reindexed descending by total cells on both axes. proportions:
    counts_ct row(tissue)-normalized to %, i.e. each tissue's own celltype
    breakdown sums to 100%. proportions_col: counts_ct column(celltype)-
    normalized to %, i.e. each celltype's own tissue breakdown sums to 100%."""
    celltype_counts = df[celltype_col].value_counts().sort_values(ascending=True)
    tissue_order_asc = df[tissue_col].value_counts(ascending=True).index.tolist()

    counts_ct = df.groupby([tissue_col, celltype_col]).size().unstack(fill_value=0)
    counts_ct = counts_ct.reindex(index=tissue_order_asc[::-1])
    counts_ct = counts_ct.reindex(columns=celltype_counts.index[::-1], fill_value=0)

    proportions = counts_ct.div(counts_ct.sum(axis=1), axis=0) * 100
    proportions_col = counts_ct.div(counts_ct.sum(axis=0), axis=1) * 100
    return celltype_counts, counts_ct, proportions, proportions_col


def plot_celltype_count_bar(celltype_counts: pd.Series, color_dict: dict, outpath: Path) -> None:
    bar_colors = [color_dict.get(ct, "#999999") for ct in celltype_counts.index]
    log_counts = np.log10(celltype_counts.values)

    plt.figure(figsize=(8, 12))
    bars = plt.barh(celltype_counts.index, log_counts, color=bar_colors)
    plt.xlabel("Number of cells (log10)")
    plt.ylabel("Cell type")
    plt.title("Number of cells per cell type")
    plt.grid(False)
    for bar, count in zip(bars, celltype_counts.values):
        x = bar.get_width()
        y = bar.get_y() + bar.get_height() / 2
        plt.text(x + 0.02, y, f"{int(count)}", va="center", ha="left", fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, bbox_inches="tight")
    plt.close()


def plot_composition_heatmaps(proportions: pd.DataFrame, proportions_col: pd.DataFrame, cmap, outdir: Path) -> None:
    plt.figure(figsize=(10, 15))
    sns.heatmap(proportions.T, cmap=cmap, linewidths=0.5, linecolor="gray",
                cbar_kws={"label": "Percent of cells"}, vmin=0, vmax=100)
    plt.ylabel("Cell type")
    plt.xlabel("Tissue")
    plt.grid(False)
    plt.title("Cell type composition per tissue")
    plt.tight_layout()
    plt.savefig(outdir / "composition_celltype_per_tissue_heatmap.pdf", bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(10, 15))
    sns.heatmap(proportions_col, cmap=cmap, linewidths=0.5, linecolor="gray",
                cbar_kws={"label": "Percent of cells"}, vmin=0, vmax=100)
    plt.ylabel("Tissue")
    plt.xlabel("Cell type")
    plt.title("Tissue composition per cell type")
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(outdir / "composition_tissue_per_celltype_heatmap.pdf", bbox_inches="tight")
    plt.close()


def compute_ranked_metric(proportions: pd.DataFrame, reducer, value_name: str, color_dict: dict) -> pd.DataFrame:
    """Reduce `proportions` (tissue-normalized %, celltype columns) to one
    value per celltype via `reducer`, rank descending, and attach color."""
    values = reducer(proportions).sort_values(ascending=False)
    out = values.reset_index()
    out.columns = ["celltype", value_name]
    out["rank"] = range(1, len(out) + 1)
    out["color"] = out["celltype"].map(color_dict).fillna("#d3d3d3")
    return out


def plot_ranked_scatter(
    df_metric: pd.DataFrame, value_col: str, ylabel: str, title: str,
    ylim: tuple, outpath: Path, labeled: bool = False, save_png: bool = False,
) -> None:
    figsize = (12, 20) if labeled and value_col == "median_pct" else (12, 8) if labeled else (6, 5)
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(
        df_metric["rank"], df_metric[value_col],
        s=35 if labeled else 20, c=df_metric["color"].tolist(),
        edgecolor="black" if labeled else None, linewidth=0.3 if labeled else 0,
    )
    if labeled:
        for _, row in df_metric.iterrows():
            ax.annotate(
                row["celltype"], xy=(row["rank"], row[value_col]), xytext=(5, 0),
                textcoords="offset points", ha="left", va="center", fontsize=7,
            )
        ax.set_xlim(0.5, len(df_metric) + 8)
    ax.set_ylim(*ylim)
    ax.set_xlabel("Cell type rank (1 = highest)")
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    if save_png:
        fig.savefig(outpath.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cell-table", required=True, help="per-cell table, tab-separated, with --tissue-col/--celltype-col columns")
    ap.add_argument("--tissue-col", default="tissue")
    ap.add_argument("--celltype-col", default="celltype")
    ap.add_argument("--celltype-colors", help="optional JSON mapping celltype name -> hex color (default: auto-assigned)")
    ap.add_argument(
        "--heatmap-colors",
        default=",".join(DEFAULT_HEATMAP_COLORS),
        help="comma-separated hex colors (low to high) for the composition heatmaps' colormap, "
        f"default: the manuscript's 'dark citrus' scale ({','.join(DEFAULT_HEATMAP_COLORS)})",
    )
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.cell_table, sep="\t")
    for col in (args.tissue_col, args.celltype_col):
        if col not in df.columns:
            raise SystemExit(f"column {col!r} not found in --cell-table")
    df = df[[args.tissue_col, args.celltype_col]].dropna().astype(str)
    print(f"{len(df)} cells x {df[args.celltype_col].nunique()} celltypes x {df[args.tissue_col].nunique()} tissues")

    color_dict = (
        json.loads(Path(args.celltype_colors).read_text())
        if args.celltype_colors else default_color_dict(df[args.celltype_col].unique())
    )

    celltype_counts, counts_ct, proportions, proportions_col = compute_counts_and_proportions(
        df, args.tissue_col, args.celltype_col
    )

    plot_celltype_count_bar(celltype_counts, color_dict, outdir / "celltype_cell_counts.pdf")

    heatmap_cmap = LinearSegmentedColormap.from_list(
        "composition_cmap", args.heatmap_colors.split(",")[::-1], N=256
    )
    plot_composition_heatmaps(proportions, proportions_col, heatmap_cmap, outdir)

    df_max = compute_ranked_metric(proportions, lambda p: p.max(axis=0), "max_pct", color_dict)
    df_third_max = compute_ranked_metric(
        proportions, lambda p: p.apply(lambda s: s.nlargest(3).iloc[-1]), "third_max_pct", color_dict
    )
    df_median = compute_ranked_metric(proportions, lambda p: p.median(axis=0), "median_pct", color_dict)
    df_gini = compute_ranked_metric(proportions, lambda p: p.apply(gini_index), "gini", color_dict)

    plot_ranked_scatter(df_max, "max_pct", "Max % of cells across tissues", "", (-10, 100),
                         outdir / "celltype_max_percent_ranked.pdf")
    plot_ranked_scatter(df_third_max, "third_max_pct", "3rd Max % of cells across tissues", "", (-10, 100),
                         outdir / "celltype_third_max_percent_ranked.pdf")
    plot_ranked_scatter(df_median, "median_pct", "Median % of cells across tissues", "", (0, 15),
                         outdir / "celltype_median_percent_ranked.pdf")
    plot_ranked_scatter(df_median, "median_pct", "Median % of cells across tissues",
                         "Cell-type median percentage across tissues", (0, 15),
                         outdir / "celltype_median_percent_ranked_labeled.pdf", labeled=True, save_png=True)
    plot_ranked_scatter(df_gini, "gini", "Gini index", "", (0, 1),
                         outdir / "celltype_gini_ranked.pdf")
    plot_ranked_scatter(df_gini, "gini", "Gini index",
                         "Tissue distribution inequality by cell type\nHigher Gini = more concentrated across tissues",
                         (0, 1), outdir / "celltype_gini_ranked_labeled.pdf", labeled=True, save_png=True)

    n_cells_df = celltype_counts.rename("n_cells").reset_index()
    n_cells_df.columns = ["celltype", "n_cells"]
    metrics = (
        df_max[["celltype", "max_pct"]]
        .merge(df_third_max[["celltype", "third_max_pct"]], on="celltype")
        .merge(df_median[["celltype", "median_pct"]], on="celltype")
        .merge(df_gini[["celltype", "gini"]], on="celltype")
        .merge(n_cells_df, on="celltype")
    )
    metrics.to_csv(outdir / "celltype_composition_metrics.tsv", sep="\t", index=False)

    print(f"done -> {outdir}")


if __name__ == "__main__":
    main()
