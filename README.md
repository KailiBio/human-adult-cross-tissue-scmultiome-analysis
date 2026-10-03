# Human Adult Cross-Tissue scMultiome Analysis

Analysis scripts for the manuscript describing a cross-tissue single-nucleus
multi-omic (paired RNA + chromatin accessibility) atlas of the adult human
body.

## Abstract

Diverse human cell types establish specialized functions through lineage- and context-specific regulatory programs. Interpreting non-coding genetic risk requires integrated multi-omic reference maps that directly connect regulatory DNA to cellular expression across human tissues. Here we present a single-nucleus multi-omic atlas comprising 459,856 transcriptomic and chromatin accessibility profiles from 21 adult human tissues and four donors, including paired measurements from 160,688 nuclei. The atlas resolves nine cell lineages, 61 broad cell types and 313 subclusters, and identifies 1,085,062 candidate cis-regulatory elements (cCREs), including 161,270 novel elements absent from ENCODE. Regulatory activity was dominated by cell identity but refined by tissue context. Joint profiling enabled 871,177 cCRE-gene associations and revealed lineage-specific regulatory architectures. Cross-tissue accessibility further identified lineage-restricted and constitutively inaccessible chromatin domains, the latter showing preferential hypomethylation across human cancers. Furthermore, we leverage this dataset to train sequence-to-function models to predict chromatin-accessibility effects for 548,656 fine-mapped variants, identifying 18,133 high-effect variants, including 1,120 broadly active variants. Models trained for eight endothelial subtypes further resolve predicted variant effects across vascular beds. Together, this atlas provides a comprehensive cellular and computational framework for interpreting regulatory sequence, context-dependent gene control, and complex trait genetics across the human body.

## Relationship to the companion pipeline repository

This repository holds the **downstream analysis** behind the manuscript's
figures and results. Data preprocessing — per-sample QC, filtering, cell
annotation, and RNA+ATAC integration into the atlas — lives in a separate,
paired repository:

- [Multiomic-workflow](https://github.com/KailiBio/Multiomic-workflow) —
  QC, cell annotation, and multiome integration pipeline (data
  preprocessing).
- **human-adult-cross-tissue-scmultiome-analysis** (this repo) — downstream
  analysis scripts for the manuscript.

## Contents

This repository is being populated with the analyses underlying the
manuscript, organized by topic as scripts are added:

- [`00_atlas_overview/`](00_atlas_overview) — Tissue→lineage→celltype→subcluster Sankey diagram; celltype×tissue composition heatmaps and ranked specificity metrics
- [`01_ccre_identification/`](01_ccre_identification) — Final cCRE calling, TSS-proximity grouping, CG/GC sequence features
- [`02_ccre_analysis/`](02_ccre_analysis) — cCRE ubiquity, marker cCREs per tissue-celltype, tissue- vs. celltype-dominant cCREs, TE enrichment
- [`03_conservation/`](03_conservation) — Alignment "triangle" matrix and phyloP/conservation aggregation profiles
- [`04_TF_regulatory/`](04_TF_regulatory) — ChromVAR motif deviation scores and TF expression vs. ChromVAR correlation
- [`05_cCRE-gene_links/`](05_cCRE-gene_links) — DORC (domains of regulatory chromatin) calling per cell type, at multiple cutoffs
- [`06_chromatin_domains/`](06_chromatin_domains) — Genome-wide chromatin domain bins, AER/ADR group definitions, cancer DNA methylation, ADR subgroup circle plot
- [`07_variant_effect_models/`](07_variant_effect_models) — seq2PRINT base model training, subtype LoRA fine-tuning, predicted variant-effect distributions, GWAS PIP threshold/enrichment testing
- [`common/`](common) — shared utility code used across topic folders (filename sanitization, default color palettes)

Each topic folder contains its own `scripts/`. Set up the
base analysis environment with `conda env create -f environment.yml`; a few
scripts additionally depend on
[scPrinter](https://github.com/buenrostrolab/scPrinter) (its own GPU stack,
installed separately) -- each topic folder's README flags which ones.

## Citation

If you use this code, please cite:

> Fan et al. *Single-Nucleus Multi-Omic Atlas Maps Regulatory Architecture and
> Non-Coding Variant Effects across Adult Human Tissues*.
> bioRxiv https://doi.org/10.64898/2026.09.25.754561

## Contact

For questions or issues, please
[open an issue](https://github.com/KailiBio/human-adult-cross-tissue-scmultiome-analysis/issues)
on GitHub.
