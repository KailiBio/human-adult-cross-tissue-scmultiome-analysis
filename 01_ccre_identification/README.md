# cCRE identification

Scripts for calling candidate cis-regulatory elements (cCREs) from
cross-tissue chromatin accessibility data and grouping them by distance to
GENCODE TSS. See [`05_cCRE-gene_links/`](../05_cCRE-gene_links) for linking
cCREs to target genes.

- `scripts/` — command-line / batch scripts
- `notebooks/` — exploratory and figure-generation notebooks

## Calling the final cCRE set

`scripts/call_ccres.py` takes a fragment-level AnnData (obs must have a
tissue-celltype grouping column and `n_fragment`) and produces the final
consensus cCRE set:

1. Calls raw MACS3 peaks per tissue-celltype cluster (`snapatac2.tl.macs3`).
2. Filters out clusters without enough sequencing depth (too few cells or
   fragments), caps each remaining cluster to its top-N peaks by q-value,
   and merges them into consensus peak regions (`snapatac2.tl.merge_peaks`).
3. Intersects consensus peaks against the raw per-cluster peaks (via
   `bedtools intersect`), keeps the best (max) q-value supporting each
   consensus peak, and thresholds on q-value for the final cCRE set.

```bash
python scripts/call_ccres.py \
  --input-h5ad /path/to/fragments.h5ad \
  --blacklist /path/to/hg38_blacklist.bed \
  --chrom-sizes /path/to/hg38.chrom.sizes \
  --outdir /path/to/output_dir
```

Outputs: `cluster_depth_and_peak_counts.tsv`, `consensus_peaks.best_qvalue.tsv`,
and the final `final_ccres.bed` / `final_ccres.tsv`. Defaults for depth and
q-value thresholds (`--min-cells 300`, `--min-fragments 1e6`,
`--top-n-peaks 180000`, `--qvalue-threshold 4`) follow those used in the
manuscript; override via CLI flags as needed.

## TSS-proximity grouping

`scripts/annotate_tss_proximity.sh` groups the final cCREs into
TSS-overlap / TSS-proximal / TSS-distal by distance to the nearest GENCODE
TSS (`bedtools closest -d`):

- distance ≤ `-d` (default 200 bp) → `TSS-overlap`
- distance ≤ `-p` (default 2,000 bp) → `TSS-proximal`
- otherwise → `TSS-distal`

```bash
scripts/annotate_tss_proximity.sh \
  -a final_ccres.bed \
  -b /path/to/gencode_tss.bed \
  -o final_ccres.tss_annotated.tsv
```

`-b` should be a BED file with one row per GENCODE TSS (e.g. extracted
from a GENCODE GTF as the 5' end of each transcript).

## CG/GC sequence features

`scripts/check_cg_content.sh` checks the CG/GC sequence features of a
given cCRE list: extracts each cCRE's sequence (`bedtools getfasta`),
computes its GC content and CpG observed/expected ratio (CpG O/E =
`(CpG_count / length) / ((C_freq + G_freq) / 2)^2`, a standard measure of
CpG depletion/enrichment relative to base composition), and optionally
joins per-cCRE mono-CG and di-CG signal from pre-built genome-wide
CpG-density bigWig tracks.

```bash
scripts/check_cg_content.sh \
  -a final_ccres.bed4 \
  -f hg38.fa \
  -o final_ccres.cg_feature.txt \
  --monocg-bw hg38_monoCG.bw \
  --dicg-bw hg38_diCG.bw
```

`--monocg-bw`/`--dicg-bw` are optional; building those genome-wide tracks
(scanning the genome FASTA for mono-/di-CG density, `bedGraphToBigWig`) is
a separate, genome-level prep step and out of scope here. Output columns:
`id`, `gc_content`, `cpg_oe`, `num_cg`, `num_c`, `num_g`, `length`, and
(if given) `mono_cg`, `di_cg`.
