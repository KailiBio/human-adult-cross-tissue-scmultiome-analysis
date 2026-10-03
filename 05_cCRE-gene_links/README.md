# cCRE–gene links

Scripts for linking cCREs to putative target genes from paired RNA+ATAC
measurements.

- `scripts/` — command-line / batch scripts
- `notebooks/` — exploratory and figure-generation notebooks

## DORC (domains of regulatory chromatin) per cell type

`scripts/run_dorc_per_celltype.py` runs DORC analysis independently per
cell type: for each cell type, subsets the matching cells from paired
RNA/ATAC AnnData and correlates peak accessibility against gene expression
within a TSS-centered window
([scPrinter](https://github.com/buenrostrolab/scPrinter)'s
`scp.dorc.fast_gene_peak_corr`), keeps significant peak-gene links, and
calls DORC genes via the J-plot (genes with an outsized number of
significant peak links).

```bash
python scripts/run_dorc_per_celltype.py \
  --rna-h5ad CZI_shared_cell.RNA.h5ad \
  --atac-peak-h5ad CZI_shared_cell.ATAC_peak.h5ad \
  --celltype-key-rna CellAnnotation_L1 \
  --celltype-key-atac celltype \
  --genome hg38 \
  --outdir /path/to/output_dir
```

`--rna-h5ad`/`--atac-peak-h5ad` must already share the same cell barcodes
(obs_names) — that pairing is a separate, upstream data-prep step.
`--genome` is `hg38` or `mm10` (drives both `scp.genome` and the matching
FigR TSS range dataset). By default both positive and negative
correlations are kept; pass `--pos-only` to restrict to positive links.

Outputs: per-celltype `results/dorc_results.{celltype}.csv` (significant
peak-gene links) and `results/dorc_genelist.{celltype}.csv` (called DORC
genes), `figures/dorc_j_plot.{celltype}.png`, and the merged
`atac_peak.withDORC.h5ad` (all per-celltype results attached under
`.uns['dorc']`).

## Re-calling DORC genes at a different cutoff

`scripts/call_dorc_genes.py` factors the J-plot gene-calling step out of
`run_dorc_per_celltype.py`, so DORC genes can be re-called at a different
`--cutoff`/`--label-top` from already-computed peak-gene correlation
tables, without re-running the (expensive) correlation step.

```bash
python scripts/call_dorc_genes.py \
  --input results/ \
  --pattern "dorc_results.*.csv" \
  --cutoff 7 \
  --outdir /path/to/output_dir
```

`--input` is a single correlation table or a directory of them (e.g.
`run_dorc_per_celltype.py`'s `results/` dir); each file's label is derived
by stripping `--label-prefix` (default `dorc_results.`) from its filename.
Outputs: per-table `results/dorc_genelist.{label}.csv` and
`figures/dorc_j_plot.{label}.png` (same as above), plus
`dorc_gene_calling_summary.txt` (label, n_peak_gene_links, n_dorc_genes).
