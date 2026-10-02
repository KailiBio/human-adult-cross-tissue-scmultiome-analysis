# cCRE analysis

Scripts for characterizing the consensus cCRE set: activity tiering,
cross-cluster ubiquity, and comparisons against external references
(ENCODE cCREs, sciATAC peaks, VISTA enhancers).

- `scripts/` — command-line / batch scripts
- `notebooks/` — exploratory and figure-generation notebooks

## cCRE ubiquity

`scripts/compute_ccre_ubiquity.py` takes the final cCRE set and the
per-cluster raw peaks it was built from (both produced by
[`01_ccre_identification/scripts/call_ccres.py`](../01_ccre_identification/scripts/call_ccres.py))
and computes, for each cCRE, the number of tissue-celltype clusters it is
active in, plus a binary cCRE x tissue-celltype activity matrix.

```bash
python scripts/compute_ccre_ubiquity.py \
  --final-ccre-bed final_ccres.bed \
  --raw-peaks-bed raw_peaks.capped.bed \
  --outdir /path/to/output_dir
```

Outputs: `ccre_ubiquity_counts.tsv` (`Peaks`, `n_active_celltypes`) and
`ccre_celltype_binary_matrix.tsv` (cCREs x tissue-celltype clusters, 0/1).
