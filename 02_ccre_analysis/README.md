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

## Marker cCREs per tissue-celltype

`scripts/identify_marker_ccres.py` takes a peak x cell AnnData (e.g. the
normalized snapATAC2 peak matrix) and, for each group in a tissue-celltype
obs column, ranks cCREs by `mean(group) - mean(rest)` and keeps the top-N as
that group's marker cCREs.

```bash
python scripts/identify_marker_ccres.py \
  --input-h5ad peak_mat.nor.h5ad \
  --groupby celltype_glue \
  --n-top 500 \
  --outdir /path/to/output_dir
```

Outputs: `marker_ccres_per_group.tsv` (group, rank, `Peaks`) and
`marker_ccres_union.txt` (union of marker cCREs across all groups).

## Tissue- vs celltype-dominant cCREs

`scripts/classify_ccre_dominance.py` takes a consensus peak-call metatable
(cCRE x tissue-celltype group, 0/1) and a group -> tissue/celltype mapping,
and classifies each cCRE by its recurrence rate: how often it's called
across a celltype's own tissue instances vs. across a tissue's own
celltypes. A cCRE is celltype-dominant if its best celltype recurrence is
`>= --tau` and at least its best tissue recurrence (ties go to celltype),
tissue-dominant if the reverse, otherwise `interaction`. Only celltypes seen
in more than `--min-tissues` tissues are used on the celltype axis; tissues
with a single celltype are excluded from the tissue axis; cCREs called in
every included group are dropped as uninformative.

```bash
python scripts/classify_ccre_dominance.py \
  --peak-metatable consensus_peak.metatable.txt \
  --group-metadata celltype_glue_group_to_tissue_celltype.tsv \
  --min-tissues 5 \
  --tau 0.5 \
  --outdir /path/to/output_dir
```

`--group-metadata` must already reflect any dataset-specific group merging
(e.g. subtype or tissue-sub-region consolidation) applied upstream. Outputs:
`ccre_dominance_classification.tsv` (`Peaks`, `category`,
`max_recur_celltype`, `max_recur_tissue`) and
`ccre_dominance_fraction_by_celltype.tsv` (per selected celltype's
celltype/tissue/interaction fraction breakdown and `n_peaks`).

## TE enrichment

`scripts/check_te_enrichment.sh` takes a cCRE info table and a
RepeatMasker-derived TE annotation BED, and computes what fraction of cCREs
overlap a TE (by each cCRE's midpoint, so a wide peak spanning several TE
copies isn't double-counted) overall, by TE class (e.g. SINE, LINE, LTR),
and by TE family within SINE/LTR/LINE (e.g. SINE/Alu, LINE/L1, plus an
"other" bucket per class), stratified by two grouping columns (e.g. a
tier/category and an ENCODE-overlap flag) and compared against each TE's
background genome-wide coverage.

```bash
scripts/check_te_enrichment.sh \
  -a CZI_CRE_final.info.txt \
  -t hg38.rmsk.group.bed \
  -g hg38.chrom.sizes \
  -o /path/to/output_dir
```

`-a`/`-t` are tab-separated with no header. `-i`/`-s`/`-e` (default columns
4/7/9) point at the cCRE id and the two stratifying columns in `-a`;
`-f`/`-c` (default columns 7/10) point at the TE family/class columns in
`-t`, e.g. `LINE/L1` and `LINE`. `--classes`/`--sine-families`/
`--ltr-families`/`--line-families` default to the classes and named
families used in the manuscript, with everything else per class bucketed
as "other" — override to test a different set. Ambiguous RepeatMasker
calls (`?`-suffixed classes) are excluded from the by-class/by-family
breakdown by default.

Outputs are headerless TSVs, wide, one row per (group, encode) pair,
matching the column layout `0-7-4d_novel_peak_TE_enrichment.ipynb` already
assigns by hand (`dat.columns = [...]`) when reading these tables:
`ccre_overlap_te.txt` (group, encode, overlapTE, total, percent),
`ccre_overlap_te.by_class.txt` (+ 2 more columns per tested class:
overlap/percent), and `ccre_overlap_te.by_family.txt` (+ 2 more columns per
named family and "other" bucket, in the order SINE[Alu, MIR, other] then
LTR[ERVK, ERV1, ERVL, ERVL-MaLR, other] then LINE[L1, L2, other]). A
`te_background_coverage.log` records each tested class/family's genome-wide
bp and percent coverage.
