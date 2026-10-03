#!/usr/bin/env python
"""Identify marker cCREs per tissue-celltype group.

For each group in `--groupby` (e.g. a combined tissue+celltype label such as
`celltype_glue`), ranks cCREs by (mean accessibility in group) - (mean
accessibility in the rest) and keeps the top `--n-top` as that group's marker
cCREs. Reports the per-group ranked lists and the union of marker cCREs
across all groups.

Consolidated from 0-7-4a2_peak_dynamic_analysis.ipynb section 1 (1-2 to 1-5).
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


def rank_group_markers(
    X, var_names, group_indices: dict, n_top: int
) -> tuple[dict, set]:
    """Rank each group's cCREs by mean(group) - mean(rest) and keep the top n_top.

    Returns (group_top_peaks, union_peaks): group_top_peaks maps group name to
    its top-N peak ids ordered by descending delta; union_peaks is the set of
    peak ids selected by at least one group.
    """
    n_cells_total = X.shape[0]
    sum_total = np.asarray(X.sum(axis=0)).ravel() if sp.issparse(X) else X.sum(axis=0)

    group_top_peaks = {}
    union_peaks = set()

    for group, idx in group_indices.items():
        n_group = len(idx)
        n_other = n_cells_total - n_group
        if n_other == 0:
            continue

        X_group = X[idx, :]
        sum_group = (
            np.asarray(X_group.sum(axis=0)).ravel() if sp.issparse(X_group) else X_group.sum(axis=0)
        )
        sum_other = sum_total - sum_group

        delta = sum_group / n_group - sum_other / n_other

        n_select = min(n_top, delta.size)
        top_idx = np.argpartition(-delta, n_select - 1)[:n_select]
        top_idx = top_idx[np.argsort(-delta[top_idx])]

        peaks = list(var_names[top_idx])
        group_top_peaks[group] = peaks
        union_peaks.update(peaks)

    return group_top_peaks, union_peaks


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-h5ad", required=True, help="peak x cell AnnData (e.g. snapATAC2 peak matrix)")
    ap.add_argument("--groupby", required=True, help="obs column giving the tissue-celltype grouping, e.g. 'celltype_glue'")
    ap.add_argument("--n-top", type=int, default=500, help="number of top marker cCREs to keep per group (default: 500)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(args.input_h5ad)
    group_indices = adata.obs.groupby(args.groupby, observed=True).indices

    print(f"ranking marker cCREs for {len(group_indices)} '{args.groupby}' groups")
    group_top_peaks, union_peaks = rank_group_markers(adata.X, adata.var_names, group_indices, args.n_top)

    long_df = pd.DataFrame(
        [
            (group, rank, peak)
            for group, peaks in group_top_peaks.items()
            for rank, peak in enumerate(peaks, start=1)
        ],
        columns=[args.groupby, "rank", "Peaks"],
    )
    long_df.to_csv(outdir / "marker_ccres_per_group.tsv", sep="\t", index=False)

    pd.Series(sorted(union_peaks), name="Peaks").to_csv(
        outdir / "marker_ccres_union.txt", index=False, header=False
    )

    print(f"  {len(union_peaks)} marker cCREs in the union across {len(group_indices)} groups")


if __name__ == "__main__":
    main()
