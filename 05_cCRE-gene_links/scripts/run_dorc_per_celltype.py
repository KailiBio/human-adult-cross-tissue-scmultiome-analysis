#!/usr/bin/env python
"""Run DORC (gene-peak correlation) analysis independently per cell type.

For each cell type in the ATAC cell-by-peak AnnData, subsets the matching
cells from a paired RNA AnnData, correlates peak accessibility against
gene expression within a TSS-centered window (scPrinter's
`scp.dorc.fast_gene_peak_corr`), keeps significant peak-gene links, and
calls "domains of regulatory chromatin" (DORC) genes via the J-plot
(`scp.dorc.dorc_j_plot`, genes with an outsized number of significant
peak links). Per-celltype results are written out individually and also
collected into the ATAC AnnData's `.uns['dorc']`.

RNA and ATAC AnnData must already share the same cell barcodes (obs_names)
-- that pairing is a separate, upstream data-prep step.

Consolidated from 1-1-1_run_DORC_perCelltype.ipynb section 1-1.
"""

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc
import scprinter as scp
from matplotlib import pyplot as plt


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rna-h5ad", required=True, help="RNA AnnData, cells x genes, aligned to --atac-peak-h5ad (same obs_names)")
    ap.add_argument("--atac-peak-h5ad", required=True, help="ATAC cell x peak AnnData, aligned to --rna-h5ad (same obs_names)")
    ap.add_argument("--celltype-key-rna", default="CellAnnotation_L1", help="obs column in --rna-h5ad with cell type labels (default: CellAnnotation_L1)")
    ap.add_argument("--celltype-key-atac", default="celltype", help="obs column in --atac-peak-h5ad with cell type labels (default: celltype)")
    ap.add_argument("--genome", default="hg38", choices=["hg38", "mm10"], help="genome build: drives both scp.genome and the matching FigR TSS range dataset (default: hg38)")
    ap.add_argument("--window-pad-size", type=int, default=50000, help="bp padding around each gene's TSS to search for correlated peaks (default: 50000)")
    ap.add_argument("--n-jobs", type=int, default=32, help="parallel jobs (default: 32)")
    ap.add_argument("--n-bg", type=int, default=100, help="number of background peaks per test, for the null distribution (default: 100)")
    ap.add_argument("--pval-cutoff", type=float, default=0.05, help="keep peak-gene links with pvalZ <= this (default: 0.05)")
    ap.add_argument("--pos-only", action="store_true", help="keep only positively correlated peak-gene links (default: both positive and negative)")
    ap.add_argument("--multimapping", action="store_true", help="allow a peak to link to multiple genes (default: off)")
    ap.add_argument("--j-plot-cutoff", type=float, default=7, help="J-plot inflection cutoff for calling a gene a DORC (default: 7)")
    ap.add_argument("--label-top", type=int, default=25, help="number of top DORC genes to label on the J-plot (default: 25)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    results_dir = outdir / "results"
    figures_dir = outdir / "figures"
    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    genome = getattr(scp.genome, args.genome)
    tss_df = getattr(scp.datasets, f"FigR_{args.genome}TSSRanges")

    adata_rna = sc.read_h5ad(args.rna_h5ad)
    adata_atac_peak = sc.read_h5ad(args.atac_peak_h5ad)

    celltypes = adata_atac_peak.obs[args.celltype_key_atac].unique().tolist()
    print(f"running DORC for {len(celltypes)} cell types")

    adata_atac_peak.uns.setdefault("dorc", {})
    for celltype in celltypes:
        print(f"processing {celltype}")
        celltype_std = celltype.replace(" ", "_").replace("/", "_")

        adata_rna_sub = adata_rna[adata_rna.obs[args.celltype_key_rna] == celltype, :].copy()
        adata_atac_sub = adata_atac_peak[adata_rna_sub.obs_names, :].copy()

        dorc_all = scp.dorc.fast_gene_peak_corr(
            adata_atac_sub,
            adata_rna_sub,
            genome=genome,
            tss_df=tss_df,
            gene_list=None,
            window_pad_size=args.window_pad_size,
            n_jobs=args.n_jobs,
            n_bg=args.n_bg,
            pval_cut=None,
            pos_only=args.pos_only,
            multimapping=args.multimapping,
        )
        dorc = dorc_all[dorc_all["pvalZ"] <= args.pval_cutoff]
        dorc.to_csv(results_dir / f"dorc_results.{celltype_std}.csv", sep="\t", index=False)

        gene_list = scp.dorc.dorc_j_plot(dorc, cutoff=args.j_plot_cutoff, label_top=args.label_top, return_gene_list=True)
        pd.Series(gene_list, name="gene").to_csv(results_dir / f"dorc_genelist.{celltype_std}.csv", index=False)
        plt.savefig(figures_dir / f"dorc_j_plot.{celltype_std}.png", dpi=300, bbox_inches="tight")
        plt.close()

        adata_atac_peak.uns["dorc"][celltype] = dorc.copy()

    adata_atac_peak.write(outdir / "atac_peak.withDORC.h5ad")
    print(f"done: {len(adata_atac_peak.uns['dorc'])} cell types")


if __name__ == "__main__":
    main()
