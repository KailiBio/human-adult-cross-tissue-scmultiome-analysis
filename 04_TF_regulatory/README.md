# Transcription factor analysis

Scripts for transcription factor motif and footprinting analysis within
cCREs.

- `scripts/` — command-line / batch scripts

## ChromVAR motif deviation scores

`scripts/run_chromvar.py` builds a cell x peak matrix from a pre-built
[scPrinter](https://github.com/buenrostrolab/scPrinter) printer object and a
peak set (e.g. the final cCRE set), then runs ChromVAR: samples GC-matched
background peaks, scans TF motifs, computes per-cell motif deviation
scores, and bags motifs with highly correlated deviations into
non-redundant TF groups.

```bash
python scripts/run_chromvar.py \
  --printer CZI_snATAC_base_scprinter.h5ad \
  --peak-bed final_ccres.bed9 \
  --genome hg38 \
  --species human \
  --device cuda \
  --outdir /path/to/output_dir
```

The printer object must already be built from fragments
(`scp.pp.import_fragments`/`scp.pp.call_peaks`) — that's a separate,
heavier data-prep step and out of scope here. `--peak-bed` is BED9
(chrom, start, end, name, score, strand, source, ...).

Outputs: `cell_peak.h5ad` (cell x peak matrix), `chromvar.h5ad` (per-cell,
per-motif deviation scores), `chromvarBG.h5ad` (deviations subset to the
bagged, non-redundant TF list), and `chromvar_bagging_tflist.txt` (the
bagging result mapping retained TFs to their correlated group).

## TF expression vs. ChromVAR correlation

`scripts/compute_tf_chromvar_correlation.py` correlates each TF's RNA
expression against its own ChromVAR motif deviation score: for every
shared TF and cell type, computes the Spearman r/R^2/p across cells of
that type (with per-celltype downsampling), then BH-FDR corrects the
p-values two ways — within each celltype (across TFs) and globally (across
every TF x celltype instance).

```bash
python scripts/compute_tf_chromvar_correlation.py \
  --rna-h5ad rna_shared.h5ad \
  --chromvar-h5ad chromvarBG.h5ad \
  --celltype-key CellAnnotation_L1 \
  --outdir /path/to/output_dir
```

`--chromvar-h5ad` should be `run_chromvar.py`'s `chromvarBG.h5ad` (the
bagged, non-redundant TF list) rather than the raw `chromvar.h5ad` --
correlating redundant, highly-correlated motifs separately just inflates
the multiple-testing burden without adding information. `--rna-h5ad`/
`--chromvar-h5ad` must already share the same cell barcodes and the same
TF/gene set — reconciling RNA cell ids with ATAC/ChromVAR barcodes
(sample/channel/donor string-matching) is a separate, cohort-specific
data-prep step done upstream of this script.

Outputs: `tf_chromvar_correlation.per_celltype.txt` (raw, pre-FDR),
`tf_chromvar_correlation.per_celltype.ctFDR.txt` (FDR within each
celltype — this is the table the volcano plots below read),
`tf_chromvar_correlation.per_celltype.globalFDR.txt` (FDR across every
instance), and `per_tf/TF_R_per_celltype.{TF}.txt` (one table per TF).

`scripts/plot_tf_chromvar_volcano.py` takes the `...ctFDR.txt` table and
makes two kinds of volcano plot (r vs. -log10(FDR)): one per TF across
cell types, and one per cell type across all TFs.

```bash
python scripts/plot_tf_chromvar_volcano.py \
  --corr-table tf_chromvar_correlation.per_celltype.ctFDR.txt \
  --celltype-colors celltype_colors.tsv \
  --outdir /path/to/output_dir
```

`--celltype-colors` is optional (tab-separated `celltype`, `color`, no
header) and only affects the per-TF volcano plots' point colors.
Per-TF plots additionally filter to `n_cells >= --min-cells-per-tf-volcano`
(default 100). Outputs: `TF_volcano/TF_volcano_by_celltype.{TF}.volcano.png|pdf`
and `celltype_volcano/volcano_allTFs.celltype_{celltype}.png|pdf`.
