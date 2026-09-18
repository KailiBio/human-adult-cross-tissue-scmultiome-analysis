# Human Adult Cross-Tissue scMultiome Analysis

Analysis scripts for the manuscript describing a cross-tissue single-nucleus
multi-omic (paired RNA + chromatin accessibility) atlas of the adult human
body.

## Abstract

A body-wide view of human gene regulation requires measuring transcription
and chromatin accessibility together in the same nuclei across tissues.
Here we generated a cross-tissue single-nucleus multi-omic atlas of 459,856
nuclei from 21 adult human tissues and four donors, with paired RNA and
chromatin profiles from the same nuclei. This scale, tissue breadth and
matched design resolve nine major lineages and 61 broad cell types,
revealing conserved cell-type programs reshaped by tissue context. The
atlas identifies 1,085,062 candidate cis-regulatory elements, including
161,270 elements absent from the ENCODE catalogue, and 871,177 cCRE–gene
links connecting chromatin accessibility to gene expression. Its
cross-tissue scope further enables a genome-wide, cell-type-resolved map of
active and repressive chromatin domains, distinguishing lineage-specific
release of repression from constitutively repressive regions shared across
tissues. These constitutive domains show preferential and tumour-type-specific
remodelling across human cancers. Finally, sequence models trained across
atlas-defined cell types identify 18,133 predicted chromatin-accessibility
high-effect variants, including 1,120 broadly active variants with
increased trait pleiotropy, and predict differential variant effects across
vascular endothelial subtypes. Together, these data establish a cross-tissue
RNA–chromatin framework linking cell identity and tissue context to
regulatory elements, target genes and disease-associated non-coding
variation.

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
manuscript, organized by topic as scripts are added, including:

- cCRE identification and cCRE–gene linking
- Cross-tissue, cell-type-resolved chromatin domain (active/repressive) mapping
- Chromatin domain remodelling analysis across human cancers
- Sequence-based models for predicting chromatin-accessibility variant effects

## Citation

If you use this code, please cite: *[citation placeholder — update once the
manuscript is published]*.

## Contact

For questions or issues, please
[open an issue](https://github.com/KailiBio/human-adult-cross-tissue-scmultiome-analysis/issues)
on GitHub.
