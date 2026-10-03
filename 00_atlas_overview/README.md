# Atlas overview

Figures summarizing the overall composition of the atlas (tissue, lineage,
celltype, and subcluster structure), independent of any specific downstream
analysis.

- `scripts/` — command-line / batch scripts
- `notebooks/` — exploratory and figure-generation notebooks

## Tissue -> lineage -> celltype -> subcluster-count Sankey diagram

`scripts/plot_tissue_lineage_celltype_sankey.py` draws a 4-column Sankey
diagram from a per-cell table: tissue -> lineage -> celltype as destination-
colored bezier ribbons between uniform-width node bars (with the tissue and
lineage stacks compressed relative to, and top-aligned with, the full-height
celltype stack), plus a 4th column of uniform-height boxes, one per celltype,
each labeled with that celltype's total number of subclusters (e.g. Leiden
clusters), summed across every tissue it appears in.

Tissue and lineage node order is chosen by a few rounds of flow-weighted
barycenter reordering (the standard two-layer crossing-reduction heuristic)
to untangle the many-to-many tissue<->lineage flows; celltype is grouped by
parent lineage (in that reordered sequence) and sorted by size within each
lineage.

```bash
python scripts/plot_tissue_lineage_celltype_sankey.py \
  --cell-table cells.tsv \
  --lineage-order "Epithelial Cell,Endothelial Cell,Mesenchymal Cell,Muscle Cell,Immune Cell,Neuroendocrine Cell,Nerve Cell,Glial Cell,Pigment Cell" \
  --outprefix figures/atlas_sankey
```

`--cell-table` is a per-cell table, tab-separated, with `--tissue-col`/
`--lineage-col`/`--celltype-col`/`--subcluster-col` columns (default:
`tissue`, `CellAnnotation_L0`, `CellAnnotation_L1`, `leiden_new` — e.g.
exported from an AnnData's `.obs`). Any desired label collapsing/splitting
(e.g. merging rare subtypes into one label, renaming multi-site tissues)
should already be applied upstream. `--lineage-order`/`--tissue-order` are
optional comma-separated preferred orders (default: alphabetical); the
manuscript's lineage order is given in the example above.
`--lineage-colors`/`--celltype-colors`/`--tissue-colors` are optional JSON
color maps (default: an auto-assigned deterministic palette).
`--multi-word-organ-prefixes` (default `Small_Intestine`) keeps tissue names
that themselves contain `_` whole when grouping tissues into organs for
inter-organ spacing. Note the subcluster count is **not** deduplicated
across tissues: a subcluster id shared by a celltype in two tissues is
counted once per tissue.

Several layout constants (node width, flow opacity, gap sizes, label
spacing/offsets, barycenter rounds) are exposed as flags, defaulting to the
values used in the manuscript figure — see `--help`.

Outputs: `<--outprefix>.pdf` and `<--outprefix>.png`.
