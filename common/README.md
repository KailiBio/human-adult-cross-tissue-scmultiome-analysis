# Common

Shared utility code used by scripts across multiple topic folders. Each
script that uses one of these adds the repo root to `sys.path` at import
time (so it works regardless of the caller's working directory) via:

```python
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common.naming import sanitize_label
```

- `naming.py` — `sanitize_label(label)`: makes a free-text label (e.g. a
  celltype or subtype name) filesystem-safe by collapsing any run of
  non-alphanumeric characters to a single `_`. Used wherever a celltype/
  subtype name is turned into part of an output filename, so the same name
  always sanitizes to the same string across scripts:
  `07_variant_effect_models/scripts/train_seq2print_per_celltype.py`,
  `07_variant_effect_models/scripts/lora_finetune_vec_subtypes.py`,
  `05_cCRE-gene_links/scripts/run_dorc_per_celltype.py`,
  `04_TF_regulatory/scripts/plot_tf_chromvar_volcano.py`.
- `colors.py` — `default_color_dict(categories)`: deterministic fallback
  color palette (tab20 + tab20b + tab20c, 60 colors, cycled) for a set of
  categories with no explicit color map. Used by
  `00_atlas_overview/scripts/plot_tissue_lineage_celltype_sankey.py` and
  `00_atlas_overview/scripts/plot_celltype_tissue_composition.py` when their
  `--*-colors` JSON argument isn't given.
