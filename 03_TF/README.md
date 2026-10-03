# Transcription factor analysis

Scripts for transcription factor motif and footprinting analysis within
cCREs.

- `scripts/` — command-line / batch scripts
- `notebooks/` — exploratory and figure-generation notebooks

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
