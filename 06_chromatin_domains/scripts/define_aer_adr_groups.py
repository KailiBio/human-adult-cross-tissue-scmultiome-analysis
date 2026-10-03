#!/usr/bin/env python
"""Define AER/ADR groups from per-bin ADR ubiquity.

Classifies each genomic bin (e.g. the 100kb bins from
get_chromatin_domain_bins.sh) into one of four ATAC-depleted-region (ADR)
groups by how many of the --n-groups cell types/clusters it was called an
ADR in ("ubiquity"):

  - AER          ubiquity == 0              (never an ADR -- accessible in every group)
  - cts-ADR      0 < ubiquity <= cts_max     (cell-type-specific ADR)
  - middle-ADR   cts_max < ubiquity < ubi_min
  - ubi-ADR      ubiquity >= ubi_min         (ubiquitous ADR, depleted almost everywhere)

cts_max/ubi_min default to 10%/90% of --n-groups (the thresholds used in
the manuscript, e.g. 16/144 of 160 cell types), rounded to the nearest
integer; override with --cts-max-pct/--ubi-min-pct.

Consolidated from 1-3-3c_ADR_activity.claude.ipynb section 2.1.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def define_groups(ubiquity: pd.Series, n_groups: int, cts_max_pct: float, ubi_min_pct: float) -> pd.Series:
    cts_max = round(cts_max_pct * n_groups)
    ubi_min = round(ubi_min_pct * n_groups)
    assert cts_max < ubi_min, "--cts-max-pct must be lower than --ubi-min-pct"

    group = pd.Series(np.full(len(ubiquity), "middle-ADR", dtype=object), index=ubiquity.index)
    group[ubiquity == 0] = "AER"
    group[(ubiquity > 0) & (ubiquity <= cts_max)] = "cts-ADR"
    group[ubiquity >= ubi_min] = "ubi-ADR"
    return group


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bed", required=True, help="per-bin table, tab-separated, no header (e.g. genome_100kb.clean.ubiquity_mapability.txt)")
    ap.add_argument("--id-col", type=int, default=4, help="1-based column with the bin id (default: 4)")
    ap.add_argument("--ubiquity-col", type=int, default=5, help="1-based column with the ADR ubiquity count (default: 5)")
    ap.add_argument("--n-groups", type=int, required=True, help="total number of cell types/clusters the ubiquity count was computed over")
    ap.add_argument("--cts-max-pct", type=float, default=0.10, help="upper ubiquity fraction (of --n-groups) for cts-ADR (default: 0.10)")
    ap.add_argument("--ubi-min-pct", type=float, default=0.90, help="lower ubiquity fraction (of --n-groups) for ubi-ADR (default: 0.90)")
    ap.add_argument("--outfile", required=True)
    args = ap.parse_args()

    outfile = Path(args.outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)

    bed = pd.read_csv(args.bed, sep="\t", header=None)
    ids = bed[args.id_col - 1]
    ubiquity = pd.to_numeric(bed[args.ubiquity_col - 1], errors="coerce").fillna(0)

    group = define_groups(ubiquity, args.n_groups, args.cts_max_pct, args.ubi_min_pct)

    out = pd.DataFrame({"id": ids, "ubiquity": ubiquity.values, "group": group.values})
    out.to_csv(outfile, sep="\t", index=False)

    counts = out["group"].value_counts().reindex(["AER", "cts-ADR", "middle-ADR", "ubi-ADR"]).fillna(0).astype(int)
    print(counts.to_string())


if __name__ == "__main__":
    main()
