# Atlas overview

Figures summarizing the overall composition of the atlas (tissue, lineage,
celltype, and subcluster structure), independent of any specific downstream
analysis.

- `scripts/` — command-line / batch scripts

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

## Celltype x tissue composition: counts, heatmaps, and ranked specificity metrics

`scripts/plot_celltype_tissue_composition.py` takes a per-cell table of
(tissue, celltype) labels and produces a cells-per-celltype bar chart, two
composition heatmaps (each tissue's celltype breakdown, and each celltype's
tissue breakdown), and four ranked scatter plots of each celltype's
distribution across tissues computed from the tissue-normalized composition
table: its single highest per-tissue percentage ("max %"), its third-highest
("3rd max %"), its median percentage, and the Gini index of that percentage
across tissues (0 = evenly spread, 1 = concentrated in one tissue). Median
and Gini each additionally get a labeled ("every point annotated with its
celltype name") version.

```bash
python scripts/plot_celltype_tissue_composition.py \
  --cell-table cells.tsv \
  --outdir /path/to/output_dir
```

`--cell-table` is tab-separated with `--tissue-col`/`--celltype-col` columns
(default: `tissue`, `celltype`). `--celltype-colors` is an optional JSON
color map (default: an auto-assigned deterministic palette).
`--heatmap-colors` is an optional comma-separated low-to-high hex color
scale for the two heatmaps (default: the manuscript's "dark citrus" scale).

Outputs: `celltype_cell_counts.pdf`, `composition_celltype_per_tissue_heatmap.pdf`,
`composition_tissue_per_celltype_heatmap.pdf`, `celltype_{max,third_max}_percent_ranked.pdf`,
`celltype_{median,gini}_ranked.pdf` plus `..._labeled.pdf`/`.png`, and
`celltype_composition_metrics.tsv` (celltype, max/3rd-max/median %, Gini,
n_cells — not written out by the source notebook, which only plotted these;
added here since the values are already computed and a table chains into
downstream scripts more easily than a figure).

Note: the source notebook's "tissue composition per cell type" heatmap had
a copy-paste bug — it plotted the transposed table under axis labels and an
inline comment that both describe the un-transposed orientation. This script
uses the un-transposed orientation (tissue rows x celltype columns),
matching the labels/comment rather than the stray `.T`.
