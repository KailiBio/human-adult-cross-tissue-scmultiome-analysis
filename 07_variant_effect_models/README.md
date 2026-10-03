# Sequence-based models for predicting chromatin-accessibility variant effects

Sequence-to-function model training and inference code used to predict
chromatin-accessibility effects for fine-mapped variants, including the
endothelial-subtype models resolving predicted effects across vascular
beds.

- `scripts/` — training / inference / evaluation scripts
- `notebooks/` — exploratory and figure-generation notebooks

## seq2PRINT base model training, per cell type

`scripts/train_seq2print_per_celltype.py` generates a
[seq2PRINT](https://github.com/buenrostrolab/scPrinter) training config
per cell type — using all of that cell type's cells as one group and its
own per-celltype peak set — and prints (or, with `--launch`, submits) the
training command for each.

```bash
python scripts/train_seq2print_per_celltype.py \
  --barcode-csv CZI_ATAC_barcode_info.BigCelltype.csv \
  --celltype-key celltype \
  --work-dir /path/to/seq2print_project \
  --genome hg38 \
  --gpus 3
```

Requires a scPrinter printer object and a seq2PRINT-preset peak set
already built per cell type under `--work-dir/celltype_printer/` and
`--work-dir/celltype_peak/` (per-celltype fragment extraction,
`scp.pp.import_fragments`, `scp.pp.call_peaks`) — a separate, heavier
data-prep step and out of scope here; cell types missing either file are
skipped with a warning. `--launch` actually submits each (expensive,
long-running GPU) training job — by default the script only prints the
command to run. Outputs: one config JSON per cell type under
`--work-dir/configs/`.

## LoRA fine-tuning a base model on its subtypes

`scripts/lora_finetune_vec_subtypes.py` builds a multi-group LoRA
fine-tuning config on top of an already-trained per-celltype base model
(one embedding-conditioned group per subtype, e.g. endothelial subtypes
within the Blood Vascular Endothelial Cell base model), then prints (or,
with `--launch`, submits) the training command.

```bash
python scripts/lora_finetune_vec_subtypes.py \
  --printer celltype_printer/CZI_snATAC_base_scprinter.Blood_Vascular_Endothelial_Cell.h5ad \
  --region-path celltype_peak/seq2print_cleaned_narrowPeak.Blood_Vascular_Endothelial_Cell.bed \
  --pretrain-model model/CZI_base.Blood_Vascular_Endothelial_Cell_fold0.pt \
  --base-model-config configs/CZI_base.Blood_Vascular_Endothelial_Cell.fold0.JSON \
  --subtype-barcodes vec_subtype_barcodes.tsv \
  --embeddings vec_rna_pca_embeddings.tsv \
  --work-dir /path/to/seq2print_project \
  --model-name CZI_LoRA.VEC \
  --gpus 1
```

`--subtype-barcodes` (columns: `barcode`, `subtype`) and `--embeddings`
(barcode-indexed, e.g. RNA PCA components) must already be prepared:
reconciling RNA subtype annotations with ATAC barcodes and filtering
subtypes by sequencing depth is cohort-specific data prep and out of
scope here. `--pretrain-model`/`--base-model-config` come from
`train_seq2print_per_celltype.py`'s output for the matching base cell
type. Despite the name (this method was originally used to fine-tune on
Blood Vascular Endothelial Cell subtypes), it's general to any base
cell type / subtype grouping. Output: one LoRA config JSON per fold under
`--work-dir/configs/`.

## Distribution of predicted variant effects

`scripts/plot_variant_effect_distribution.py` characterizes the
distribution of a seq2PRINT model's predicted chromatin-accessibility
variant effects: the pooled and per-cell-type distribution of delta
effects, each variant's cell-type specificity (Gini index of |delta|
across cell types — 0 = spread evenly, 1 = concentrated in one cell
type), and the joint density of specificity vs. effect magnitude.

```bash
python scripts/plot_variant_effect_distribution.py \
  --snp-delta SNP_delta.all.tsv \
  --celltype-colors celltype_colors.json \
  --outdir /path/to/output_dir
```

`--snp-delta` is tab-separated with a `variant_id` column and one
`*_delta` column per cell type (a `V2.` prefix and the `_delta` suffix are
stripped to recover the cell type name). `--celltype-colors` is optional
(JSON mapping cell type name to color, for the per-celltype violin plots;
default: a single color for all). Outputs: four figures
(`distribution_delta_effect_pooled`, `..._per_celltype`, `..._gini`,
`gini_vs_max_delta`, each `.png`/`.pdf`) and
`variant_effect_specificity.tsv` (`variant_id`, `delta_gini`,
`n_nonzero_celltypes`, `max_abs_delta`).
