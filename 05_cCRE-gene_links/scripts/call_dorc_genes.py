#!/usr/bin/env python
"""Call DORC genes from already-computed peak-gene correlation results.

Takes one or more DORC peak-gene correlation tables (already filtered to
significant links, e.g. the per-celltype `dorc_results.*.csv` files from
run_dorc_per_celltype.py) and calls "domains of regulatory chromatin"
(DORC) genes via the J-plot (scPrinter's `scp.dorc.dorc_j_plot`): a gene is
called a DORC if it has >= `--cutoff` significant peak links. Factored out
from the (expensive) peak-gene correlation computation so genes can be
re-called at a different cutoff without recomputing correlations.

Consolidated from 1-1-1_run_DORC_perCelltype.ipynb section 1-1 (the
DORC-gene-calling step).
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scprinter as scp


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, help="a single DORC results table, or a directory of them (see --pattern)")
    ap.add_argument("--pattern", default="dorc_results.*.csv", help="glob pattern to match DORC result files when --input is a directory (default: dorc_results.*.csv)")
    ap.add_argument("--label-prefix", default="dorc_results.", help="prefix stripped from each matched filename's stem to derive its label (default: dorc_results.)")
    ap.add_argument("--sep", default="\t", help="field separator of the input table(s) (default: tab, matching run_dorc_per_celltype.py's output)")
    ap.add_argument("--cutoff", type=float, default=7, help="a gene is called DORC if it has >= this many significant peak links (default: 7)")
    ap.add_argument("--label-top", type=int, default=25, help="number of top DORC genes to label on the J-plot (default: 25)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    results_dir = outdir / "results"
    figures_dir = outdir / "figures"
    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    input_path = Path(args.input)
    files = sorted(input_path.glob(args.pattern)) if input_path.is_dir() else [input_path]
    if not files:
        raise SystemExit(f"no files matched {args.pattern!r} under {input_path}")

    print(f"calling DORC genes (cutoff={args.cutoff}) for {len(files)} table(s)")
    summary_rows = []
    for f in files:
        label = f.stem
        if label.startswith(args.label_prefix):
            label = label[len(args.label_prefix):]
        print(f"processing {label}")

        dorc = pd.read_csv(f, sep=args.sep)
        gene_list = scp.dorc.dorc_j_plot(dorc, cutoff=args.cutoff, label_top=args.label_top, return_gene_list=True)
        pd.Series(gene_list, name="gene").to_csv(results_dir / f"dorc_genelist.{label}.csv", index=False)
        plt.savefig(figures_dir / f"dorc_j_plot.{label}.png", dpi=300, bbox_inches="tight")
        plt.close()

        summary_rows.append({"label": label, "n_peak_gene_links": len(dorc), "n_dorc_genes": len(gene_list)})

    pd.DataFrame(summary_rows).to_csv(outdir / "dorc_gene_calling_summary.txt", sep="\t", index=False)
    print("done")


if __name__ == "__main__":
    main()
