# Conservation

Scripts for evaluating sequence conservation (e.g. phyloP) across cCREs
and chromatin domains.

- `scripts/` — command-line / batch scripts
- `notebooks/` — exploratory and figure-generation notebooks

## Alignment "triangle" matrix

`scripts/generate_alignment_triangle_matrix.sh` takes per-species genome
alignment bedGraphs (e.g. a Zoonomia-style multi-species alignment) and a
cCRE set, and for each species computes the per-cCRE mean alignment score
(via `bigWigAverageOverBed`). Per cCRE, it then tallies how many species
have alignment ≥ `--align-high` ("align90") vs. ≤ `--align-low`
("align10") — the x/y axes of the "triangle plot" (named for the
triangular point cloud `align90 + align10 <= n_species` produces), used to
flag accelerated/non-conserved accessible elements.

```bash
scripts/generate_alignment_triangle_matrix.sh \
  -a final_ccres.bed4 \
  -b zoonomia_alignment/bg \
  -g hg38.chrom.sizes \
  -o /path/to/output_dir \
  -r random_windows.bed4 \
  -i final_ccres.info.txt \
  -p hg38_PAR.bed
```

`-b` is a directory of one `<species>.bg` alignment bedGraph per species;
bigWigs are built once into `-b/bw` (or `--bw-dir`) and reused on later
runs. `-r` optionally scores a background/random-region set the same way,
for comparison. `-i` (optional) is a BED-like cCRE info table
(chrom/start/end in the first 3 columns, so it can be filtered by `-p`)
used to stratify the cCRE triangle table by group/group2 (default columns
4/5/6 for id/group/group2 — override with `--id-col`/`--group-col`/
`--group2-col`); without `-i`, the cCRE table is just (id, align90,
align10), matching the background table's shape.

Outputs: `alignment_score_matrix.{cre,background}.txt` (id x species mean
alignment score) and `alignment_triangle_counts.{cre,background}.txt` (the
triangle-plot-ready align90/align10 tallies).

## phyloP/conservation aggregation profile

`scripts/generate_phylop_aggregation.sh` generates a conservation
aggregation profile for a given cCRE list: recenters each cCRE to a
fixed-width window around its midpoint, splits that window into equal
bins, computes the mean bigWig signal (phyloP, phastCons, GERP, ...) per
bin per cCRE (`bigWigAverageOverBed`), then reduces across all cCREs to a
single mean profile ready for plotting.

```bash
scripts/generate_phylop_aggregation.sh \
  -a final_ccres.bed4 \
  -w hg38_240mammal_phyloP.bw \
  -o phyloP241.CRE.mean.txt \
  -n 400 \
  --window 2000 \
  --label CRE
```

It's single-set and general by design: to compare against a background or
a stratified subset (e.g. TSS-distal only, ENCODE-overlapping only), run
it again on a random-region BED or a pre-filtered cCRE BED with a
different `--label`, then combine the resulting profile rows yourself
(e.g. `cat profile1.txt profile2.txt > combined.txt`) for plotting.
`--chrom-sizes` additionally clips windows that run past a chromosome's
end (the original script only clipped the start). cCRE ids must not
contain `_` (ids get suffixed `_<bin index>` internally and split back
apart — this is checked up front and errors clearly if violated).

Outputs: the mean profile at `-o` (one row: label, then one mean-signal
column per bin), and, with `--keep-matrix`, the per-cCRE x bin signal
matrix at `<-o>.matrix.txt`.
