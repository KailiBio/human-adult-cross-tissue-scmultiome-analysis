#!/usr/bin/env python
"""4-column Sankey diagram: tissue -> lineage -> celltype -> per-celltype
subcluster-count box.

Columns 1-3 (tissue -> lineage -> celltype) are drawn as destination-colored
bezier ribbons between uniform-width node bars, with the tissue and lineage
node stacks compressed to a fraction of the (full-height) celltype stack and
top-aligned to it -- both reduce how much vertical space the few/large
tissue and lineage categories take up relative to the many/small celltype
categories. Tissue and lineage node order is chosen by a few rounds of
flow-weighted barycenter reordering (the standard two-layer crossing-
reduction heuristic) to untangle the many-to-many tissue<->lineage flows;
celltype is simply grouped by parent lineage (in that reordered sequence)
and sorted by size within each lineage. Column 4 adds one uniform-height box
per celltype (not part of the actual flow -- an annotation only), labeled
with the number of distinct subclusters (e.g. Leiden clusters) found under
that celltype, summed across every tissue it appears in: a subcluster id
shared by one celltype across several tissues is counted once per tissue,
not deduplicated to a single cross-tissue total.

Consolidated from 0-9-3_make_final_figs.ipynb section 3g (reusing section
3's barycenter tissue/lineage ordering, 3b's compressed/top-aligned
tissue+lineage layout, and 3e's destination-colored ribbon/label style for
columns 1-3).
"""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch


def organ_group_fn(multi_word_prefixes):
    """Collapse a tissue name to its parent organ (text before the first '_'),
    except for prefixes in `multi_word_prefixes` that themselves contain '_'
    and must be kept whole (e.g. "Small_Intestine_Duodenum" -> "Small_Intestine",
    not "Small")."""

    def organ_group(tissue):
        for p in multi_word_prefixes:
            if tissue.startswith(p):
                return p
        return tissue.split("_")[0]

    return organ_group


def barycenter(target_cats, neighbor_rank, flow_series, target_is_left):
    """Mean rank (in `neighbor_rank`) of each target category's connected
    partners, weighted by flow size -- the standard two-layer crossing-
    reduction heuristic."""
    bary = {}
    for cat in target_cats:
        try:
            sub = flow_series.xs(cat, level=0 if target_is_left else 1)
        except KeyError:
            bary[cat] = np.inf
            continue
        ranks = np.array([neighbor_rank[c] for c in sub.index])
        bary[cat] = float(np.average(ranks, weights=sub.values))
    return bary


def stack_labels(ys, min_sep):
    """Nudge a sorted list of target y's apart so labels don't overlap, then
    settle back toward the true positions from both ends (simple two-pass
    solver)."""
    ys = list(ys)
    for i in range(1, len(ys)):
        if ys[i] - ys[i - 1] < min_sep:
            ys[i] = ys[i - 1] + min_sep
    for i in range(len(ys) - 2, -1, -1):
        if ys[i + 1] - ys[i] < min_sep:
            ys[i] = ys[i + 1] - min_sep
    return ys


def stacked_positions(order, sizes, gap_between):
    """Stack `order` top-to-bottom in raw cell-count units; `gap_between(prev,
    cat)` gives the gap (in cell-count units) to insert before each category
    after the first, so groups can get bigger gaps than within-group
    neighbors."""
    y = 0.0
    pos = {}
    prev = None
    for cat in order:
        if prev is not None:
            y += gap_between(prev, cat)
        h = sizes[cat]
        pos[cat] = (y, y + h)
        y += h
        prev = cat
    return pos


def compress_and_align(pos_raw, raw_max, ref_max, height_frac, align="top"):
    """Rescale a column's raw (cell-count-unit) positions so its total stack
    occupies only `height_frac` of `ref_max` -- proportions between its own
    nodes are preserved. `align` controls where the leftover slack goes:
    "top" (flush with the reference column, slack at the bottom) or "center".
    Also returns the scale factor, needed to size ribbons leaving/entering
    this column consistently with its compressed nodes."""
    scale = height_frac * ref_max / raw_max
    offset = ref_max * (1 - height_frac) / 2 if align == "center" else 0.0
    positions = {cat: (offset + y0 * scale, offset + y1 * scale) for cat, (y0, y1) in pos_raw.items()}
    return positions, scale


def default_color_dict(categories):
    """Deterministic fallback palette (tab20+tab20b+tab20c, 60 colors, cycled)
    when no explicit color map is given for a column."""
    colors = []
    for name in ("tab20", "tab20b", "tab20c"):
        colors.extend(mcolors.rgb2hex(c) for c in plt.get_cmap(name).colors)
    return {cat: colors[i % len(colors)] for i, cat in enumerate(sorted(categories))}


def order_categories(preferred_order, present):
    """`preferred_order` filtered to categories actually present, with any
    present-but-unlisted categories appended (alphabetically); falls back to
    plain alphabetical order if no preference is given."""
    if preferred_order:
        order = [c for c in preferred_order if c in present]
        missing = sorted(present - set(order))
        return order + missing
    return sorted(present)


def bezier_ribbon(x0, x1, y0_l, y1_l, y0_r, y1_r, buffer_frac=0.15):
    """Flat stubs at each end (so the ribbon visibly hugs its node) connected
    by a two-control-point bezier curve through the x midpoint."""
    buf = (x1 - x0) * buffer_frac
    xa, xb = x0 + buf, x1 - buf
    xm = (xa + xb) / 2
    verts = [
        (x0, y0_l), (xa, y0_l),
        (xm, y0_l), (xm, y0_r), (xb, y0_r),
        (x1, y0_r),
        (x1, y1_r), (xb, y1_r),
        (xm, y1_r), (xm, y1_l), (xa, y1_l),
        (x0, y1_l), (x0, y0_l),
    ]
    codes = [
        MplPath.MOVETO, MplPath.LINETO,
        MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
        MplPath.LINETO,
        MplPath.LINETO, MplPath.LINETO,
        MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
        MplPath.LINETO, MplPath.CLOSEPOLY,
    ]
    return MplPath(verts, codes)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--cell-table",
        required=True,
        help="per-cell table, tab-separated, with --tissue-col/--lineage-col/--celltype-col/--subcluster-col "
        "columns (e.g. exported from an AnnData's .obs); any desired label collapsing/splitting (e.g. merging "
        "rare subtypes, renaming multi-site tissues) should already be applied upstream",
    )
    ap.add_argument("--tissue-col", default="tissue")
    ap.add_argument("--lineage-col", default="CellAnnotation_L0")
    ap.add_argument("--celltype-col", default="CellAnnotation_L1")
    ap.add_argument(
        "--subcluster-col",
        default="leiden_new",
        help="fine-grained subcluster id column (e.g. Leiden cluster); column 4 labels each celltype with the "
        "number of distinct values of this column, summed per tissue -- not deduplicated across tissues (default: leiden_new)",
    )
    ap.add_argument("--lineage-colors", help="optional JSON mapping lineage name -> hex color (default: auto-assigned)")
    ap.add_argument("--celltype-colors", help="optional JSON mapping celltype name -> hex color (default: auto-assigned)")
    ap.add_argument("--tissue-colors", help="optional JSON mapping tissue name -> hex color (default: auto-assigned)")
    ap.add_argument(
        "--lineage-order",
        help="optional comma-separated preferred top-to-bottom lineage order (default: alphabetical). The "
        "manuscript's order was 'Epithelial Cell,Endothelial Cell,Mesenchymal Cell,Muscle Cell,Immune Cell,"
        "Neuroendocrine Cell,Nerve Cell,Glial Cell,Pigment Cell'",
    )
    ap.add_argument("--tissue-order", help="optional comma-separated preferred top-to-bottom tissue order (default: alphabetical)")
    ap.add_argument(
        "--multi-word-organ-prefixes",
        default="Small_Intestine",
        help="comma-separated tissue-name prefixes that themselves contain '_' and should be kept whole when "
        "grouping tissues into organs for the between-organ gap spacing (default: Small_Intestine, e.g. so "
        "Small_Intestine_Duodenum isn't split into organ 'Small')",
    )
    ap.add_argument("--barycenter-rounds", type=int, default=6, help="rounds of flow-weighted barycenter reordering for tissue/lineage (default: 6)")
    ap.add_argument("--gap-frac", type=float, default=0.0015, help="within-group node gap, as a fraction of total cells (default: 0.0015)")
    ap.add_argument("--between-group-gap-mult", type=float, default=8.0, help="multiple of the within-group gap used between tissue organs / between every lineage (default: 8)")
    ap.add_argument("--compressed-height-frac", type=float, default=0.40, help="tissue's and lineage's total stack height, as a fraction of the (uncompressed) celltype column's span (default: 0.40)")
    ap.add_argument("--node-width", type=float, default=0.045, help="node/box width, shared by all 4 columns (default: 0.045)")
    ap.add_argument("--flow-alpha", type=float, default=0.65, help="flow ribbon opacity (default: 0.65)")
    ap.add_argument("--label-gap", type=float, default=0.05, help="near label offset from its node (default: 0.05)")
    ap.add_argument("--label-gap-far", type=float, default=0.32, help="far (second-depth) label offset for celltype's staggered leader lines (default: 0.32)")
    ap.add_argument("--min-label-sep-frac", type=float, default=0.005, help="minimum vertical clearance between stacked labels, as a fraction of the celltype column's full span (default: 0.005)")
    ap.add_argument("--figsize", default="16,16", help="figure size in inches, WIDTH,HEIGHT (default: 16,16)")
    ap.add_argument("--outprefix", required=True, help="output path prefix; writes <prefix>.pdf and <prefix>.png")
    args = ap.parse_args()

    df = pd.read_csv(args.cell_table, sep="\t")
    required_cols = [args.tissue_col, args.lineage_col, args.celltype_col, args.subcluster_col]
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        raise SystemExit(f"column(s) not found in --cell-table: {missing_cols}")

    df = df[required_cols].dropna(subset=[args.tissue_col, args.lineage_col, args.celltype_col]).copy()
    df[args.tissue_col] = df[args.tissue_col].astype(str)
    df[args.lineage_col] = df[args.lineage_col].astype(str)
    df[args.celltype_col] = df[args.celltype_col].astype(str)
    if df[args.subcluster_col].isna().any():
        raise SystemExit(f"--cell-table has missing values in {args.subcluster_col!r} -- every cell must have a subcluster id")

    total = len(df)
    print(
        f"{total} cells x {df[args.lineage_col].nunique()} lineages x "
        f"{df[args.celltype_col].nunique()} celltypes x {df[args.tissue_col].nunique()} tissues"
    )

    organ_group = organ_group_fn(args.multi_word_organ_prefixes.split(","))

    lineage_pref = args.lineage_order.split(",") if args.lineage_order else None
    tissue_pref = args.tissue_order.split(",") if args.tissue_order else None
    lineage_order = order_categories(lineage_pref, set(df[args.lineage_col].unique()))
    tissue_order = order_categories(tissue_pref, set(df[args.tissue_col].unique()))

    l1_to_l0 = dict(df[[args.celltype_col, args.lineage_col]].drop_duplicates().values)
    tissue_counts = df[args.tissue_col].value_counts()
    flow_tissue_lineage = df.groupby([args.tissue_col, args.lineage_col]).size()

    # --- barycenter crossing-reduction: alternately reorder tissue and lineage
    # by the flow-weighted mean rank of what they connect to (tissue grouped by
    # organ throughout, so e.g. Colon_Sigmoid/_Transverse stay adjacent) ---
    for _ in range(args.barycenter_rounds):
        tissue_rank = {c: i for i, c in enumerate(tissue_order)}
        bary_lineage = barycenter(lineage_order, tissue_rank, flow_tissue_lineage, target_is_left=False)
        lineage_order = sorted(lineage_order, key=lambda c: bary_lineage[c])

        lineage_rank = {c: i for i, c in enumerate(lineage_order)}
        bary_tissue = barycenter(tissue_order, lineage_rank, flow_tissue_lineage, target_is_left=True)
        organ_members = {}
        for t in tissue_order:
            organ_members.setdefault(organ_group(t), []).append(t)
        organ_bary = {
            organ: np.average([bary_tissue[t] for t in members], weights=[tissue_counts[t] for t in members])
            for organ, members in organ_members.items()
        }
        tissue_order = sorted(tissue_order, key=lambda t: (organ_bary[organ_group(t)], bary_tissue[t]))

    # celltype: grouped by parent lineage (in the settled lineage_order), then by size within lineage
    lineage_rank = {c: i for i, c in enumerate(lineage_order)}
    celltype_counts = df[args.celltype_col].value_counts()
    celltype_order = sorted(
        df[args.celltype_col].unique(),
        key=lambda c: (lineage_rank[l1_to_l0[c]], -celltype_counts[c]),
    )

    color_dicts = {
        "tissue": json.loads(Path(args.tissue_colors).read_text()) if args.tissue_colors else default_color_dict(df[args.tissue_col].unique()),
        "lineage": json.loads(Path(args.lineage_colors).read_text()) if args.lineage_colors else default_color_dict(df[args.lineage_col].unique()),
        "celltype": json.loads(Path(args.celltype_colors).read_text()) if args.celltype_colors else default_color_dict(df[args.celltype_col].unique()),
    }

    # --- compressed+top-aligned tissue/lineage stacks, full-height celltype stack ---
    gap = total * args.gap_frac
    between_group_gap = gap * args.between_group_gap_mult

    tissue_sizes = df[args.tissue_col].value_counts()
    tissue_pos_raw = stacked_positions(
        tissue_order, tissue_sizes,
        lambda prev, cat: gap if organ_group(prev) == organ_group(cat) else between_group_gap,
    )
    tissue_raw_max = max(y1 for y0, y1 in tissue_pos_raw.values())

    lineage_sizes = df[args.lineage_col].value_counts()
    lineage_pos_raw = stacked_positions(lineage_order, lineage_sizes, lambda prev, cat: between_group_gap)
    lineage_raw_max = max(y1 for y0, y1 in lineage_pos_raw.values())

    celltype_sizes = df[args.celltype_col].value_counts()
    celltype_pos_raw = stacked_positions(celltype_order, celltype_sizes, lambda prev, cat: gap)
    celltype_raw_max = max(y1 for y0, y1 in celltype_pos_raw.values())
    ref_max = celltype_raw_max

    tissue_positions, tissue_scale = compress_and_align(tissue_pos_raw, tissue_raw_max, ref_max, args.compressed_height_frac)
    lineage_positions, lineage_scale = compress_and_align(lineage_pos_raw, lineage_raw_max, ref_max, args.compressed_height_frac)
    celltype_positions, celltype_scale = compress_and_align(celltype_pos_raw, celltype_raw_max, ref_max, 1.0)

    col_positions = [tissue_positions, lineage_positions, celltype_positions]
    col_scale = [tissue_scale, lineage_scale, celltype_scale]
    col_order = [tissue_order, lineage_order, celltype_order]
    col_key = [args.tissue_col, args.lineage_col, args.celltype_col]
    col_role = ["tissue", "lineage", "celltype"]

    node_width = args.node_width
    x_positions = np.arange(4)

    # --- columns 1-3 ribbons: destination-colored (lineage color for
    # tissue->lineage, celltype color for lineage->celltype) ---
    cursor_out = [dict((cat, pos[0]) for cat, pos in col_positions[ci].items()) for ci in range(3)]
    cursor_in = [dict((cat, pos[0]) for cat, pos in col_positions[ci].items()) for ci in range(3)]

    ribbons = []
    for ci in range(2):
        left_col, right_col = col_key[ci], col_key[ci + 1]
        flow_counts = df.groupby([left_col, right_col]).size()
        color_dict_ci = color_dicts[col_role[ci + 1]]
        for left_cat in col_order[ci]:
            for right_cat in col_order[ci + 1]:
                if (left_cat, right_cat) not in flow_counts.index:
                    continue
                n = int(flow_counts.loc[(left_cat, right_cat)])

                y0_l = cursor_out[ci][left_cat]
                y1_l = y0_l + n * col_scale[ci]
                cursor_out[ci][left_cat] = y1_l

                y0_r = cursor_in[ci + 1][right_cat]
                y1_r = y0_r + n * col_scale[ci + 1]
                cursor_in[ci + 1][right_cat] = y1_r

                x0 = x_positions[ci] + node_width / 2
                x1 = x_positions[ci + 1] - node_width / 2
                path = bezier_ribbon(x0, x1, y0_l, y1_l, y0_r, y1_r)
                color = color_dict_ci.get(right_cat, "#999999")
                ribbons.append((n, path, color))

    # --- column 4: uniform-height box per celltype (annotation-only, no real
    # flow), labeled with its total subcluster count ---
    subcluster_counts_by_celltype_tissue = df.groupby([args.celltype_col, args.tissue_col])[args.subcluster_col].nunique()
    subcluster_counts_by_celltype = subcluster_counts_by_celltype_tissue.groupby(level=0).sum()

    box_gap = gap
    n_celltypes = len(celltype_order)
    box_height = (ref_max - box_gap * (n_celltypes - 1)) / n_celltypes
    celltype_box_positions = {}
    y = 0.0
    for cat in celltype_order:
        celltype_box_positions[cat] = (y, y + box_height)
        y += box_height + box_gap

    for cat in celltype_order:
        y0_l, y1_l = celltype_positions[cat]
        y0_r, y1_r = celltype_box_positions[cat]
        x0 = x_positions[2] + node_width / 2
        x1 = x_positions[3] - node_width / 2
        path = bezier_ribbon(x0, x1, y0_l, y1_l, y0_r, y1_r)
        color = color_dicts["celltype"].get(cat, "#999999")
        ribbons.append((int(celltype_sizes[cat]), path, color))

    # --- draw ---
    figsize = tuple(float(v) for v in args.figsize.split(","))
    fig, ax = plt.subplots(figsize=figsize)

    # smallest ribbons first, largest last (on top), so large flows don't bury
    # the many smaller ones crossing through the same region
    for n, path, color in sorted(ribbons, key=lambda r: r[0]):
        ax.add_patch(PathPatch(path, facecolor=color, edgecolor="none", alpha=args.flow_alpha, lw=0))

    for ci in range(3):
        color_dict = color_dicts[col_role[ci]]
        x = x_positions[ci]
        for cat in col_order[ci]:
            y0, y1 = col_positions[ci][cat]
            color = color_dict.get(cat, "#999999")
            ax.add_patch(plt.Rectangle((x - node_width / 2, y0), node_width, y1 - y0,
                                        facecolor=color, edgecolor="black", linewidth=0.8, zorder=3))

    x4 = x_positions[3]
    for cat, (y0, y1) in celltype_box_positions.items():
        color = color_dicts["celltype"].get(cat, "#999999")
        ax.add_patch(plt.Rectangle((x4 - node_width / 2, y0), node_width, y1 - y0,
                                    facecolor=color, edgecolor="black", linewidth=0.8, zorder=3))

    min_label_sep = ref_max * args.min_label_sep_frac

    # tissue and lineage: plain beside-node labels, no leader line (few enough
    # categories that de-overlapping is the only thing needed)
    for ci in (0, 1):
        x = x_positions[ci]
        side_left = (ci == 0)
        anchor_x = x - node_width / 2 if side_left else x + node_width / 2
        label_x = anchor_x - args.label_gap if side_left else anchor_x + args.label_gap
        items = sorted(((y0 + y1) / 2, cat) for cat, (y0, y1) in col_positions[ci].items())
        stacked_ys = stack_labels([y for y, _ in items], min_label_sep)
        for (true_y, cat), label_y in zip(items, stacked_ys):
            ax.text(label_x, label_y, cat, fontsize=6, ha="right" if side_left else "left",
                    va="center", zorder=5, color="black")

    # celltype: staggered leader-line callouts (dense cluster of rare celltypes
    # needs the extra room a leader line buys)
    x = x_positions[2]
    anchor_x = x + node_width / 2
    items = sorted(((y0 + y1) / 2, cat) for cat, (y0, y1) in col_positions[2].items())
    stacked_ys = stack_labels([y for y, _ in items], min_label_sep)
    for i, ((true_y, cat), label_y) in enumerate(zip(items, stacked_ys)):
        label_x = anchor_x + (args.label_gap if i % 2 == 0 else args.label_gap_far)
        ax.plot([anchor_x, label_x], [true_y, label_y], color="#888", lw=0.5, zorder=4)
        ax.text(label_x, label_y, cat, fontsize=6, ha="left", va="center", zorder=5, color="black")

    # column 4: one label per uniform-height box, "<celltype> (<n_subclusters>)"
    x = x_positions[3]
    anchor_x = x + node_width / 2
    for cat, (y0, y1) in celltype_box_positions.items():
        n_sub = int(subcluster_counts_by_celltype.get(cat, 0))
        ax.text(anchor_x + args.label_gap, (y0 + y1) / 2, f"{cat} ({n_sub})", fontsize=6,
                ha="left", va="center", zorder=5, color="black")

    ax.set_xlim(-0.9, 3.9)
    ax.set_ylim(-ref_max * 0.05, ref_max * 1.02)
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()

    outprefix = Path(args.outprefix)
    outprefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{outprefix}.pdf", bbox_inches="tight")
    fig.savefig(f"{outprefix}.png", dpi=200, bbox_inches="tight")
    print(f"done -> {outprefix}.pdf / .png")


if __name__ == "__main__":
    main()
