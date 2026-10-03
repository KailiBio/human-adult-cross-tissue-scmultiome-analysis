# Cross-tissue, cell-type-resolved chromatin domain mapping

Scripts for mapping active/repressive chromatin domains across cell types
and tissues, including identification of lineage-restricted and
constitutively inaccessible domains.

- `scripts/` — command-line / batch scripts

## Genome bins for chromatin domain calling

`scripts/get_chromatin_domain_bins.sh` tiles the genome into fixed-size
bins (default 100kb) as the base unit for calling chromatin domains, then
filters out bins with low mappability and/or overlapping a blacklist, so
downstream domain-calling only considers "clean" genomic bins.

```bash
scripts/get_chromatin_domain_bins.sh \
  -g hg38.chrom.sizes.clean \
  -m k50.Umap.MultiTrackMappability.bw \
  -o genome_100kb.clean.bed4 \
  -w 100000 \
  -t 0.5 \
  -b hg38_blacklist.bed \
  -b encBlacklistV2.bed \
  -b grcExclusions.bed
```

`-b` is repeatable (e.g. the ENCODE blacklist, unusual regions, and GRC
exclusion tracks used in the manuscript); converting any UCSC bigBed
blacklist tracks to BED (`bigBedToBed`) is a separate, one-time reference
prep step and out of scope here — without any `-b`, only the mappability
filter is applied. Output: a BED4 of clean bins (chrom, start, end,
`chrom:start-end`).

## Defining AER/ADR groups

`scripts/define_aer_adr_groups.py` classifies each genomic bin into one of
four ATAC-depleted-region (ADR) groups by how many cell types/clusters it
was called an ADR in ("ubiquity"): `AER` (never an ADR — accessible in
every group), `cts-ADR` (cell-type-specific), `middle-ADR`, or `ubi-ADR`
(ubiquitous, depleted almost everywhere). The `cts-ADR`/`ubi-ADR` cutoffs
default to 10%/90% of the total number of groups (the thresholds used in
the manuscript, e.g. 16/144 of 160 cell types).

```bash
python scripts/define_aer_adr_groups.py \
  --bed genome_100kb.clean.ubiquity_mapability.txt \
  --n-groups 160 \
  --outfile adr_groups.txt
```

`--bed` is tab-separated with no header (e.g. the per-bin table produced
further downstream of `get_chromatin_domain_bins.sh`'s ADR-activity
tallying), with `--id-col`/`--ubiquity-col` (default 4/5) pointing at the
bin id and ubiquity count columns. Output: `id`, `ubiquity`, `group`.

## Cancer DNA methylation at AER/ADR regions

Two scripts test whether AER/ADR regions are differentially methylated in
TCGA tumor vs. normal samples. Both expect per-sample methylation
summaries already computed from raw arrays (IDATs -> beta values, joined
against a probe/region-to-ADR map) — that's a separate, heavier data-prep
step and out of scope here.

`scripts/compute_dnam_tumor_normal_delta.py` pools all probes within each
region-type/subgroup per sample, then computes Mann-Whitney U tumor-vs-normal
tests and `delta = mean(Tumor) - mean(Normal)` per TCGA project — both raw
and background-normalized (each sample's value minus its own whole-array
global mean beta, to control for bulk hypomethylation that varies by
cancer type).

```bash
python scripts/compute_dnam_tumor_normal_delta.py \
  --batch-summaries-glob "results/pilot_raw/*.out.csv" \
  --batch-summaries-glob "results/full_raw/*.out.csv" \
  --outdir /path/to/output_dir
```

`--batch-summaries-glob` is repeatable (e.g. one per processing run/cohort
batch); each CSV has one row per sample with `global_mean`,
`{region_type}_mean`/`_n` per region type, `sg_{subgroup}_mean`/`_n` per
subgroup. Region types tested default to `AER,cts-ADR,ubi-ADR`
(`--region-types`); subgroups are discovered from the data automatically.
Outputs: `all_projects_summary.csv`, per-project Mann-Whitney stats
(`all_projects_stats_*.csv`), and heatmap-ready wide tables
(`heatmap_table_*delta*.csv`).

`scripts/compute_dnam_region_level_delta.py` instead computes delta
independently for each individual ADR region, then aggregates per-region
deltas within each subgroup (mean or median across regions) — giving every
region equal weight regardless of its probe count, unlike the pooled-probe
approach above.

```bash
python scripts/compute_dnam_region_level_delta.py \
  --adr-bed ADR.sorted.bed \
  --region-detail-glob "results/*_raw_detail/*.regions.parquet" \
  --sample-summary all_projects_summary.csv \
  --outdir /path/to/output_dir
```

`--adr-bed` maps `region_id` to `region_type`/`subgroup` (chr, start, end,
region_id, celltype_desc, numeric_code, region_type, subgroup).
`--region-detail-glob` is repeatable, matching per-sample per-region mean
beta (long format: sample_uuid, region_id, mean_beta, n). Outputs:
`region_level_delta*.csv` (per-project, per-region deltas) and
`heatmap_table_subgroup_delta_by_region*.csv` (mean/median-aggregated,
raw/global-normalized).

## Genome-wide ADR subgroup circle plot

`scripts/plot_adr_subgroup_circle_plot.py` draws a circos-style, whole-genome
plot of 100kb-bin chromatin subgroup assignments: chromosome ideogram, gene
density (discrete color bins), domain accessibility (continuous heatmap --
inverted ADR ubiquity, so higher = more open), and one consolidated presence
ring per region group (which bins were assigned to any subgroup in that
group, adjacent bins merged into runs). Matches the "v4" style of the source
notebook: a chosen chromosome dropped from the ideogram entirely (default
chrY), the circle rotated so a chosen chromosome sits at the 3 o'clock
position (default chr16), and the accessibility ring drawn as a per-bin
heatmap rather than a bar chart (which rasterized incorrectly in PDF export).

```bash
python scripts/plot_adr_subgroup_circle_plot.py \
  --subgroup-bed ADR_subgroup_assignments_allRegions.bed \
  --gene-density-bed genome_100kb.geneDensity.GENCODEv47.bed4 \
  --activity-bed genome_100kb.clean.ubiquity_mapability.group.txt \
  --n-groups 160 \
  --outprefix figures/adr_circle_plot
```

`--subgroup-bed` is a per-bin table (tab-separated, no header, 8 columns:
chrom, start, end, id, label, num, group, subgroup) -- the subgroup
assignment itself (clustering cts-ADR bins by which cell type most often
calls them an ADR, e.g. `T_Cell_1`/`BVEC`/`Neuron_1`) is a separate,
upstream step and out of scope here. `--activity-bed` is the same kind of
per-bin ubiquity table [`define_aer_adr_groups.py`](#defining-aeradr-groups)
consumes (built from [`get_chromatin_domain_bins.sh`](#genome-bins-for-chromatin-domain-calling)'s
clean bins); `--n-groups`/`--activity-ubiquity-col` match that script's
flags of the same name. `--palette` (subgroup -> color) and `--region-groups`
(ring order outer-to-inner, mapping a group label to its member subgroups
and color) are optional JSON overrides of the manuscript's defaults -- see
`--help`. `--drop-chromosomes`/`--rotate-to-chrom` default to `chrY`/`chr16`
(empty string to disable either). Requires `pycirclize` and `lxml`.

Outputs: `<--outprefix>.png`, `.pdf` (ring patches rasterized at high DPI so
Illustrator gets lightweight embedded bitmaps instead of tens of thousands
of tiny per-bin paths; everything else -- labels, ticks, legend -- stays
vector), and `.svg` (every ring/label's patches grouped under one `gid` each,
so each ring is one selectable, editable object in Illustrator).
