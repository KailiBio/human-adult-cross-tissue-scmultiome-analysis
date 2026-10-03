#!/usr/bin/env python
"""Classify cCREs as tissue-dominant, celltype-dominant, or interaction.

Given a consensus peak-call metatable (cCRE x tissue-celltype group, 0/1) and
a group -> tissue/celltype mapping, computes each cCRE's recurrence rate
within every celltype (fraction of that celltype's own tissue instances the
cCRE is called in) and within every tissue (fraction of that tissue's
celltypes the cCRE is called in, using the tissue's full cellular
composition). A cCRE is:

  - tissue-dominant if its best tissue recurrence >= `--tau` and >= its best
    celltype recurrence
  - celltype-dominant if its best celltype recurrence >= `--tau` and >=
    its best tissue recurrence (ties go to celltype)
  - interaction otherwise

Only celltypes seen in more than `--min-tissues` tissues are used on the
celltype axis (not enough cross-tissue replication otherwise); tissues with
only one celltype are excluded from the tissue axis (their cCREs would be
trivially "100% tissue-recurrent"). cCREs called in every included group are
dropped as uninformative.

Assumes `--peak-metatable` and `--group-metadata` already have one column/row
per final tissue-celltype group (apply any dataset-specific group merging,
e.g. subtype or tissue-sub-region consolidation, upstream of this script).

Consolidated from 0-7-4a2_peak_dynamic_analysis.ipynb section 3.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def classify_dominance(
    peak_meta: pd.DataFrame,
    group_meta: pd.DataFrame,
    min_tissues: int,
    tau: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Classify each cCRE as celltype-, tissue-dominant, or interaction.

    peak_meta: cCRE x group boolean table. group_meta: group-indexed table
    with 'tissue' and 'celltype' columns, aligned to peak_meta's columns.
    Returns (dominant, by_celltype): per-cCRE classification (plus the
    recurrence rates behind it), and the celltype/tissue/interaction fraction
    breakdown for each well-replicated celltype's own peak set.
    """
    group_meta = group_meta.loc[peak_meta.columns]

    n_tissue_per_celltype = group_meta.groupby("celltype", observed=True)["tissue"].nunique()
    sel_celltypes = (
        n_tissue_per_celltype[n_tissue_per_celltype > min_tissues]
        .sort_values(ascending=False)
        .index.tolist()
    )
    glue_sub = group_meta[group_meta["celltype"].isin(sel_celltypes)]

    peak_meta_sub = peak_meta[glue_sub.index]
    active = peak_meta_sub.any(axis=1)
    peak_meta_sub = peak_meta_sub.loc[active].astype(np.float32)

    peak_sets_by_celltype = {
        ct: peak_meta.index[peak_meta[glue_sub.index[glue_sub["celltype"] == ct]].any(axis=1)]
        for ct in sel_celltypes
    }

    X = peak_meta_sub.values  # cCREs x groups (selected celltypes only)
    celltype_of_group = glue_sub.loc[peak_meta_sub.columns, "celltype"].astype(str)
    C = pd.get_dummies(celltype_of_group).values.astype(np.float32)
    celltype_mean = (X @ C) / C.sum(axis=0, keepdims=True)

    # tissue-side recurrence uses ALL celltypes present in a tissue, not just
    # the selected ones -- a tissue's cross-celltype consistency should
    # reflect its full cellular composition
    peak_meta_active = peak_meta.loc[peak_meta_sub.index]
    tissue_of_group_all = group_meta.loc[peak_meta_active.columns, "tissue"].astype(str)
    T_all = pd.get_dummies(tissue_of_group_all).values.astype(np.float32)
    tissue_mean = (peak_meta_active.values.astype(np.float32) @ T_all) / T_all.sum(axis=0, keepdims=True)

    full_coverage = pd.crosstab(group_meta["tissue"], group_meta["celltype"])
    n_celltypes_per_tissue = (full_coverage > 0).sum(axis=1)
    singleton_tissues = n_celltypes_per_tissue[n_celltypes_per_tissue < 2].index.tolist()
    tissue_mask = ~pd.get_dummies(tissue_of_group_all).columns.isin(singleton_tissues)

    max_recur_celltype = celltype_mean.max(axis=1)
    max_recur_tissue = tissue_mean[:, tissue_mask].max(axis=1)

    # drop cCREs called in literally every included group -- trivially uninformative
    keep = X.sum(axis=1) < X.shape[1]
    peak_index = peak_meta_sub.index[keep]
    max_recur_celltype, max_recur_tissue = max_recur_celltype[keep], max_recur_tissue[keep]

    ct_flag = max_recur_celltype >= tau
    ts_flag = max_recur_tissue >= tau
    category = np.full(len(peak_index), "interaction", dtype=object)
    category[ts_flag & (max_recur_tissue >= max_recur_celltype)] = "tissue"
    category[ct_flag & (max_recur_celltype >= max_recur_tissue)] = "celltype"  # ties -> celltype

    dominant = pd.DataFrame(
        {
            "category": category,
            "max_recur_celltype": max_recur_celltype,
            "max_recur_tissue": max_recur_tissue,
        },
        index=peak_index,
    )

    rows, n_peaks_used = {}, {}
    for ct in sel_celltypes:
        peaks_ct = peak_sets_by_celltype[ct].intersection(dominant.index)
        n_peaks_used[ct] = peaks_ct.size
        rows[ct] = (
            dominant.loc[peaks_ct, "category"]
            .value_counts(normalize=True)
            .reindex(["celltype", "tissue", "interaction"])
            .fillna(0)
        )

    by_celltype = pd.DataFrame(rows).T.reindex(sel_celltypes)
    by_celltype["n_peaks"] = pd.Series(n_peaks_used)

    return dominant, by_celltype


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--peak-metatable",
        required=True,
        help="cCRE x tissue-celltype-group consensus peak-call table (0/1), tab-separated, cCRE id in first column",
    )
    ap.add_argument(
        "--group-metadata",
        required=True,
        help="tab-separated table with group id in first column (matching --peak-metatable's columns) and 'tissue', 'celltype' columns",
    )
    ap.add_argument(
        "--min-tissues",
        type=int,
        default=5,
        help="a celltype needs >this many distinct tissues to be used on the celltype axis (default: 5)",
    )
    ap.add_argument(
        "--tau",
        type=float,
        default=0.5,
        help="recurrence-rate threshold for calling a cCRE tissue- or celltype-dominant (default: 0.5)",
    )
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    peak_meta = pd.read_csv(args.peak_metatable, sep="\t", index_col=0).astype(bool)
    group_meta = pd.read_csv(args.group_metadata, sep="\t", index_col=0)[["tissue", "celltype"]]

    print(f"classifying {peak_meta.shape[0]} cCREs across {peak_meta.shape[1]} tissue-celltype groups (tau={args.tau})")
    dominant, by_celltype = classify_dominance(peak_meta, group_meta, args.min_tissues, args.tau)

    dominant.to_csv(outdir / "ccre_dominance_classification.tsv", sep="\t", index_label="Peaks")
    by_celltype.to_csv(outdir / "ccre_dominance_fraction_by_celltype.tsv", sep="\t")

    print(dominant["category"].value_counts())


if __name__ == "__main__":
    main()
