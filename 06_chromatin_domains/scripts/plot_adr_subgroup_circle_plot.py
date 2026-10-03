#!/usr/bin/env python
"""Genome-wide circos-style plot of ADR chromatin subgroups (v4 style).

Rings, outer to inner: chromosome ideogram (cytoband), gene density
(discrete bins, sequential color ramp), domain accessibility (continuous
heatmap; higher = more open/less ADR-ubiquitous), and one consolidated
presence ring per region group (e.g. AER / Neuron cts-ADR / ubi-ADR -- which
100kb bins were assigned to any subgroup in that group, merged into runs).

Matches "v4" of the source notebook: a chosen chromosome is dropped from the
ideogram entirely (default chrY), the circle is rotated so a chosen
chromosome sits at the 3 o'clock position (default chr16), and every ring is
drawn as a per-bin colored rect (a heatmap), including domain accessibility
-- an earlier bar-ring version of that ring rasterized incorrectly in PDF
export.

Requires `pycirclize` and `lxml`.

Consolidated from 1-8_ADR_circle_plot.ipynb's "V4" section (cells 20-21).
Note: the source notebook defines a `brewer_marine` palette and documents it
as the gene-density ring's color scheme (also used for that ring's label
text, "#0868AC"), but the actual gene_cmap/legend swatches were built from a
different, never-renamed `brewer_green` variable. This script uses
brewer_marine, matching the documented intent and the label color.
"""

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lxml import etree
from pycirclize import Circos
from pycirclize import config as pycirclize_config
from pycirclize.parser import Bed
from pycirclize.utils import load_eukaryote_example_dataset

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"

SVG_NS = "http://www.w3.org/2000/svg"

DEFAULT_PALETTE = {
    "AER": "#9E9E9E",
    "T_Cell_1": "#4D9221",
    "BVEC": "#FFBE0B",
    "Glial_Cell": "#795548",
    "Macrophage": "#00bfae",
    "T_Cell_2": "#a4ff5a",
    "LEC": "#ff5b13",
    "Theca_Cell": "#e72c80",
    "Neuron_1": "#65c4ff",
    "Neuron_2": "#0e3aa9",
    "ubi-ADR": "#8e02da",
}
DEFAULT_REGION_GROUPS = {
    "AER": {"subgroups": ["AER"], "color": "#9E9E9E"},
    "Neuron cts-ADR": {"subgroups": ["Neuron_1", "Neuron_2"], "color": "#3a7fd4"},
    "ubi-ADR": {"subgroups": ["ubi-ADR"], "color": "#8e02da"},
}
BREWER_MARINE = ["#F7FCF0", "#E0F3DB", "#CCEBC5", "#A8DDB5", "#7BCCC4",
                 "#4EB3D3", "#2B8CBE", "#0868AC", "#084081"]
BREWER_HEAT = ["#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59",
               "#EF6548", "#D7301F", "#B30000", "#7F0000"]


def group_svg_by_gid(in_path, out_path, prefixes=("ring_", "label_")):
    """Collect all SVG elements sharing a gid (with one of the given
    prefixes) into a single <g> group, so each ring/label is one selectable,
    editable object in Illustrator instead of thousands of loose paths."""
    tree = etree.parse(in_path)
    root = tree.getroot()

    seen, names = set(), []
    for el in root.iter():
        gid = el.get("id")
        if gid and any(gid.startswith(p) for p in prefixes) and gid not in seen:
            seen.add(gid)
            names.append(gid)

    for name in names:
        matches = root.findall(f'.//{{{SVG_NS}}}g[@id="{name}"]')
        if not matches:
            continue
        parent = matches[0].getparent()
        group = etree.SubElement(parent, f"{{{SVG_NS}}}g")
        group.set("id", name)
        for i, el in enumerate(matches):
            el.getparent().remove(el)
            el.set("id", f"{name}__{i}")
            group.append(el)

    tree.write(out_path, xml_declaration=True, encoding="UTF-8", standalone=False)


def merge_selected(sub: pd.DataFrame) -> pd.DataFrame:
    """Merge adjacent (end[i] == start[i+1]) 100kb bins -- regardless of
    which subgroup each individual bin belongs to -- into single genomic
    intervals, per chromosome."""
    out = []
    for chrom, g in sub.sort_values("start").groupby("chrom", sort=False):
        g = g.reset_index(drop=True)
        run_id = (g["start"] != g["end"].shift()).cumsum()
        out.append(g.groupby(run_id).agg(chrom=("chrom", "first"), start=("start", "min"), end=("end", "max")))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=["chrom", "start", "end"])


def parse_r(s: str) -> tuple:
    a, b = s.split(",")
    return float(a), float(b)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--subgroup-bed", required=True,
        help="ADR subgroup assignments, tab-separated, no header, 8 columns: chrom, start, end, id, label, num, "
        "group, subgroup (100kb-bin resolution; subgroup assignment itself is a separate, upstream clustering "
        "step and out of scope here)",
    )
    ap.add_argument(
        "--gene-density-bed", required=True,
        help="gene density, tab-separated, no header, BED4: chrom, start, end, gene_count (100kb bins)",
    )
    ap.add_argument(
        "--activity-bed", required=True,
        help="per-bin ADR ubiquity table, tab-separated, no header (e.g. 06_chromatin_domains/"
        "get_chromatin_domain_bins.sh's clean bins joined against a per-bin ubiquity tally, the same table "
        "define_aer_adr_groups.py consumes)",
    )
    ap.add_argument("--activity-ubiquity-col", type=int, default=5, help="1-based column in --activity-bed with the raw ADR ubiquity count (default: 5, matching define_aer_adr_groups.py's --ubiquity-col)")
    ap.add_argument("--n-groups", type=int, required=True, help="total number of cell types/clusters the ubiquity count in --activity-bed was computed over (same meaning as define_aer_adr_groups.py's --n-groups); accessibility = --n-groups minus the raw ubiquity count, so higher = more open")
    ap.add_argument("--palette", help=f"optional JSON mapping subgroup name -> hex color (default: the manuscript's palette, {len(DEFAULT_PALETTE)} subgroups)")
    ap.add_argument(
        "--region-groups",
        help='optional JSON object, in ring order outer-to-inner, mapping group label -> {"subgroups": [...], '
        '"color": "#hex"} (default: AER / Neuron cts-ADR / ubi-ADR, matching the manuscript\'s v4 figure)',
    )
    ap.add_argument("--gene-bin-edges", default="-1,0,1,3,5,9,inf", help="comma-separated bin edges for the discrete gene-density ring (default: -1,0,1,3,5,9,inf)")
    ap.add_argument("--gene-bin-labels", default="0,1,2-3,4-5,6-9,10+", help="comma-separated labels for --gene-bin-edges (default: 0,1,2-3,4-5,6-9,10+)")
    ap.add_argument("--genome-build", default="hg38", help="genome build passed to pycirclize's bundled UCSC chrom-sizes/cytoband loader (default: hg38)")
    ap.add_argument("--chrom-bed", help="override chrom-sizes BED instead of --genome-build's bundled one")
    ap.add_argument("--cytoband-file", help="override cytoband file instead of --genome-build's bundled one")
    ap.add_argument("--drop-chromosomes", default="chrY", help="comma-separated chromosomes to drop from the ideogram entirely (default: chrY; empty string to keep all)")
    ap.add_argument("--rotate-to-chrom", default="chr16", help="rotate the circle so this chromosome's center sits at the 3 o'clock position (default: chr16; empty string to disable rotation)")
    ap.add_argument("--space", type=float, default=1.5, help="angular gap (degrees) between adjacent sectors (default: 1.5)")
    ap.add_argument("--r-gene", type=parse_r, default="88,96", help="gene-density ring (inner,outer) radius (default: 88,96)")
    ap.add_argument("--r-accessibility", type=parse_r, default="78,86", help="accessibility ring (inner,outer) radius (default: 78,86)")
    ap.add_argument("--region-ring-top", type=float, default=76, help="outer radius of the first (outermost) region-group presence ring (default: 76)")
    ap.add_argument("--region-ring-band-width", type=float, default=8, help="radial width of each region-group presence ring (default: 8)")
    ap.add_argument("--region-ring-gap", type=float, default=2, help="radial gap between consecutive region-group presence rings (default: 2)")
    ap.add_argument("--figsize", default="6,6", help="figure size in inches, WIDTH,HEIGHT (default: 6,6)")
    ap.add_argument("--png-dpi", type=int, default=300, help="PNG export DPI (default: 300)")
    ap.add_argument(
        "--pdf-raster-dpi", type=int, default=2400,
        help="DPI used to rasterize ring patches (cytoband/gene-density/accessibility/region-presence) for the "
        "PDF export only -- high because each ring draws one patch per 100kb bin, spread around the whole "
        "genome's circumference, so bins are well under 1px wide at typical DPI (default: 2400). Chromosome "
        "labels, ring-name labels, ticks, and the gene-density legend stay vector/text.",
    )
    ap.add_argument("--title", default="Genome-wide ADR regions: gene density, accessibility & consolidated subgroups", help="figure suptitle (the genome build / chromosome set is appended automatically)")
    ap.add_argument("--outprefix", required=True, help="output path prefix; writes <prefix>.png, <prefix>.pdf, <prefix>.svg")
    args = ap.parse_args()

    outprefix = Path(args.outprefix)
    outprefix.parent.mkdir(parents=True, exist_ok=True)

    palette = json.loads(Path(args.palette).read_text()) if args.palette else DEFAULT_PALETTE
    region_groups = json.loads(Path(args.region_groups).read_text()) if args.region_groups else DEFAULT_REGION_GROUPS
    gene_bin_edges = [np.inf if v == "inf" else float(v) for v in args.gene_bin_edges.split(",")]
    gene_bin_labels = args.gene_bin_labels.split(",")
    drop_chromosomes = [c for c in args.drop_chromosomes.split(",") if c]

    # ---------- load subgroup assignments, gene density, accessibility ----------
    df = pd.read_csv(
        args.subgroup_bed, sep="\t", header=None,
        names=["chrom", "start", "end", "id", "label", "num", "group", "subgroup"],
    )
    missing_colors = set(df["subgroup"].unique()) - set(palette.keys())
    if missing_colors:
        raise SystemExit(f"subgroups without a palette color: {missing_colors}")

    gene_density = pd.read_csv(args.gene_density_bed, sep="\t", header=None,
                                names=["chrom", "start", "end", "gene_count"])

    activity_raw = pd.read_csv(args.activity_bed, sep="\t", header=None)
    activity = activity_raw.iloc[:, [0, 1, 2, args.activity_ubiquity_col - 1]].copy()
    activity.columns = ["chrom", "start", "end", "ubiquity"]
    activity["activity"] = args.n_groups - activity["ubiquity"]

    df_annot = (
        df.merge(gene_density, on=["chrom", "start", "end"], how="left")
        .merge(activity[["chrom", "start", "end", "activity"]], on=["chrom", "start", "end"], how="left")
    )
    n_missing = df_annot[["gene_count", "activity"]].isna().sum().sum()
    if n_missing:
        raise SystemExit(f"{n_missing} bins in --subgroup-bed had no matching --gene-density-bed/--activity-bed row")

    activity_vmin, activity_vmax = 0, df_annot["activity"].max()
    activity_cmap = mcolors.LinearSegmentedColormap.from_list("activity_heat", BREWER_HEAT, N=256)
    activity_norm = mcolors.Normalize(vmin=activity_vmin, vmax=activity_vmax)

    gene_cmap = mcolors.LinearSegmentedColormap.from_list("gene_marine", BREWER_MARINE, N=256)
    gene_colors = [mcolors.to_hex(gene_cmap(x)) for x in np.linspace(0, 1, len(gene_bin_labels))]
    gene_color_map = dict(zip(gene_bin_labels, gene_colors))
    df_annot["gene_bin"] = pd.cut(df_annot["gene_count"], bins=gene_bin_edges, labels=gene_bin_labels)
    df_annot["gene_color"] = df_annot["gene_bin"].map(gene_color_map)

    region_sets = {
        label: merge_selected(df[df["subgroup"].isin(spec["subgroups"])])
        for label, spec in region_groups.items()
    }
    print("region group sizes:", {label: len(r) for label, r in region_sets.items()})

    # ---------- chromosome sectors: drop chromosomes, build gap, rotate ----------
    chr_bed_file, cytoband_file, _ = (
        load_eukaryote_example_dataset(args.genome_build) if not (args.chrom_bed and args.cytoband_file)
        else (args.chrom_bed, args.cytoband_file, None)
    )
    if args.chrom_bed:
        chr_bed_file = args.chrom_bed
    if args.cytoband_file:
        cytoband_file = args.cytoband_file

    chr_sizes = pd.read_csv(chr_bed_file, sep="\t", comment="#", header=None, names=["chrom", "start", "end", "name"])
    if drop_chromosomes:
        chr_sizes = chr_sizes[~chr_sizes["chrom"].isin(drop_chromosomes)]
    chr_bed_file_filtered = f"{outprefix}.chrom_sizes.filtered.bed"
    chr_sizes.to_csv(chr_bed_file_filtered, sep="\t", header=False, index=False)

    probe = Circos.initialize_from_bed(chr_bed_file_filtered, space=args.space)
    last_sector_name = probe.sectors[-1].name
    last_sector_deg_size = probe.get_sector(last_sector_name).deg_size
    gap_last_first = args.space + last_sector_deg_size / 2
    n_sectors = len(probe.sectors)
    space_list = [args.space] * (n_sectors - 1) + [gap_last_first]

    start_deg = 0.0
    if args.rotate_to_chrom:
        probe_spaced = Circos.initialize_from_bed(chr_bed_file_filtered, space=space_list)
        rotate_mid = sum(probe_spaced.get_sector(args.rotate_to_chrom).deg_lim) / 2
        rotation = (90 - rotate_mid) % 360
        start_deg = rotation - 360 if rotation > 0 else 0.0

    circos = Circos.initialize_from_bed(chr_bed_file_filtered, start=start_deg, end=start_deg + 360, space=space_list)
    cytoband_records = Bed(cytoband_file).records

    # ---------- ring radii ----------
    r_gene, r_accessibility = args.r_gene, args.r_accessibility
    region_ring_radii = {}
    for i, label in enumerate(region_groups):
        top = args.region_ring_top - i * (args.region_ring_band_width + args.region_ring_gap)
        region_ring_radii[label] = (top - args.region_ring_band_width, top)

    for sector in circos.sectors:
        chrom = sector.name
        sub = df_annot[df_annot["chrom"] == chrom].sort_values("start")

        sector.text(chrom, r=113, size=10, weight="bold", gid="label_chrom_names")

        cytoband_track = sector.add_track((97, 100))
        cytoband_track.axis()
        for rec in cytoband_records:
            if rec.chr == chrom:
                color = pycirclize_config.CYTOBAND_COLORMAP.get(str(rec.score), "white")
                cytoband_track.rect(rec.start, rec.end, fc=color, ec="none", gid="ring_cytoband", zorder=2)
        cytoband_track.xticks_by_interval(
            50_000_000, label_formatter=lambda v: f"{v / 1e6:.0f}",
            outer=True, tick_length=1, label_size=6,
        )

        gene_track = sector.add_track(r_gene)
        gene_track.axis(fc="#ffffff", ec="#dddddd", lw=0.4)
        for _, row in sub.iterrows():
            gene_track.rect(row["start"], row["end"], fc=row["gene_color"], ec="none",
                             gid="ring_gene_density", zorder=2)

        activity_track = sector.add_track(r_accessibility)
        activity_track.axis(fc="#ffffff", ec="#dddddd", lw=0.4)
        for _, row in sub.iterrows():
            activity_track.rect(row["start"], row["end"], fc=activity_cmap(activity_norm(row["activity"])),
                                 ec="none", gid="ring_domain_activity", zorder=2)

        for label, spec in region_groups.items():
            track = sector.add_track(region_ring_radii[label])
            track.axis(fc="#ffffff", ec="#dddddd", lw=0.4)
            rows = region_sets[label]
            rows = rows[rows["chrom"] == chrom]
            gid = "ring_" + label.lower().replace(" ", "_").replace("-", "_")
            for _, row in rows.iterrows():
                track.rect(row["start"], row["end"], fc=spec["color"], ec="none", gid=gid, zorder=2)

    # ---------- ring-name labels in the widened last-sector -> first-sector gap ----------
    first_sector_name = circos.sectors[0].name
    last_end = circos.get_sector(last_sector_name).deg_lim[1]
    first_start = circos.get_sector(first_sector_name).deg_lim[0]
    first_start_unwrapped = first_start + 360 if first_start < last_end else first_start
    gap_mid_deg = (last_end + first_start_unwrapped) / 2 % 360

    def gap_label(text, r, color="#555555", gid=None):
        circos.text(text, r=r, deg=gap_mid_deg, size=6.5, color=color, ha="center", va="center", gid=gid)

    gap_label("Gene density", sum(r_gene) / 2, color=BREWER_MARINE[7], gid="label_ring_gene_density")
    gap_label("Accessibility", sum(r_accessibility) / 2, color=BREWER_HEAT[-1], gid="label_ring_domain_activity")
    for label, spec in region_groups.items():
        gid = "label_ring_" + label.lower().replace(" ", "_").replace("-", "_")
        gap_label(label, sum(region_ring_radii[label]) / 2, color=spec["color"], gid=gid)

    figsize = tuple(float(v) for v in args.figsize.split(","))
    fig = circos.plotfig(figsize=figsize, dpi=args.png_dpi)

    gene_handles = [mpatches.Patch(facecolor=c, edgecolor="none", label=l) for l, c in zip(gene_bin_labels, gene_colors)]
    fig.legend(handles=gene_handles, loc="lower left", bbox_to_anchor=(0.06, 0.06),
               fontsize=7, frameon=False, title="Genes/100kb", title_fontsize=8,
               handlelength=1.2, handleheight=1.2)

    cax = fig.add_axes([0.08, 0.17, 0.018, 0.12])
    sm = plt.cm.ScalarMappable(cmap=activity_cmap, norm=activity_norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("Accessibility", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    cbar.ax.text(0.5, -0.06, "closed", transform=cbar.ax.transAxes, ha="center", va="top", fontsize=6)
    cbar.ax.text(0.5, 1.04, "open", transform=cbar.ax.transAxes, ha="center", va="bottom", fontsize=6)

    chrom_set_desc = f"{args.genome_build}, " + ("chr1-22+X" if drop_chromosomes == ["chrY"] else "custom chromosome set")
    fig.suptitle(f"{args.title} ({chrom_set_desc})", y=1.0, fontsize=14)

    fig.savefig(f"{outprefix}.png", dpi=args.png_dpi, bbox_inches="tight")

    # rasterize ring patches for the PDF export only, so Illustrator opens one
    # lightweight embedded bitmap per ring instead of tens of thousands of tiny
    # paths, while labels/ticks/legend stay real vector/text objects
    ring_patches = [p for ax in fig.axes for p in ax.patches if (gid := p.get_gid()) and gid.startswith("ring_")]
    for p in ring_patches:
        p.set_rasterized(True)
    fig.savefig(f"{outprefix}.pdf", dpi=args.pdf_raster_dpi, bbox_inches="tight")
    for p in ring_patches:
        p.set_rasterized(False)

    svg_path = f"{outprefix}.svg"
    fig.savefig(svg_path, bbox_inches="tight")
    group_svg_by_gid(svg_path, svg_path)

    print(f"done -> {outprefix}.png / .pdf / .svg")


if __name__ == "__main__":
    main()
