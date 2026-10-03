#!/usr/bin/env python
"""Run ChromVAR motif deviation scores on a cCRE x cell matrix.

Builds a cell x peak matrix from a pre-built scPrinter printer object and a
peak set, then runs ChromVAR (via scPrinter's `chromvar` module): samples
GC-matched background peaks, scans TF motifs, computes per-cell motif
deviation scores, and bags together motifs with highly correlated
deviations into non-redundant TF groups.

Requires a scPrinter printer object already built from fragments (see
`scp.pp.import_fragments`/`scp.pp.call_peaks`) -- building that object is a
separate, heavier data-prep step and out of scope here.

Consolidated from seq2print/1-1_base_model.ipynb section 1 (1-1, 1-2).
"""

import argparse
import warnings
from pathlib import Path

import scprinter as scp


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--printer", required=True, help="scPrinter printer h5ad path")
    ap.add_argument(
        "--peak-bed",
        required=True,
        help="peak set to build the cell x peak matrix from (BED9: chrom, start, end, name, "
        "score, strand, source, ...), e.g. the final cCRE set",
    )
    ap.add_argument("--region-width", type=int, default=300, help="peak width for the cell x peak matrix (default: 300)")
    ap.add_argument("--genome", default="hg38", help="scprinter.genome attribute name (default: hg38)")
    ap.add_argument("--species", default="human", choices=["human", "mouse"], help="FigR motif set to scan (default: human)")
    ap.add_argument("--device", default="cuda", choices=["cuda", "cpu"], help="device for ChromVAR deviation computation (default: cuda)")
    ap.add_argument("--gpu-device", type=int, default=0, help="GPU device id to register with rmm when --device cuda (default: 0)")
    ap.add_argument("--n-jobs", type=int, default=100, help="parallel jobs for motif scanning (default: 100)")
    ap.add_argument("--pvalue", type=float, default=5e-5, help="motif-match p-value threshold (default: 5e-5)")
    ap.add_argument("--n-bg-iterations", type=int, default=250, help="background-peak sampling iterations (default: 250)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genome = getattr(scp.genome, args.genome)
    printer = scp.load_printer(args.printer, genome)

    print(f"building cell x peak matrix from {args.peak_bed}")
    adata = scp.pp.make_peak_matrix(
        printer,
        regions=args.peak_bed,
        region_width=args.region_width,
        cell_grouping=None,
        group_names=None,
        sparse=True,
    )
    adata.write(outdir / "cell_peak.h5ad")
    printer.close()

    # only keep peaks with nonzero coverage
    coverage = adata.X.sum(axis=0)
    adata = adata[:, coverage > 0]

    if args.device == "cuda":
        warnings.filterwarnings("ignore")
        import cupy as cp
        import rmm
        from rmm.allocators.cupy import rmm_cupy_allocator

        rmm.reinitialize(managed_memory=True, pool_allocator=True, devices=args.gpu_device)
        cp.cuda.set_allocator(rmm_cupy_allocator)

    print("sampling background peaks")
    scp.chromvar.sample_bg_peaks(adata, genome=genome, method="chromvar", niterations=args.n_bg_iterations)

    print(f"scanning {args.species} motifs")
    motif_fn = scp.motifs.FigR_Human_Motifs if args.species == "human" else scp.motifs.FigR_Mouse_Motifs
    motif = motif_fn(genome, bg=list(adata.uns["bg_freq"]), n_jobs=args.n_jobs, pvalue=args.pvalue, mode="motifmatchr")
    motif.prep_scanner(None, pvalue=args.pvalue)
    motif.chromvar_scan(adata)

    print("computing per-cell motif deviations")
    chromvar = scp.chromvar.compute_deviations(adata, chunk_size=50000, device=args.device)
    chromvar.write(outdir / "chromvar.h5ad")

    print("bagging correlated motifs into TF groups")
    bagging_matrix = (
        scp.datasets.FigR_motifs_bagging_human if args.species == "human" else scp.datasets.FigR_motifs_bagging_mouse
    )
    chromvar_bg = scp.chromvar.bag_deviations(adata=chromvar, motif_motif_matrix=bagging_matrix)
    chromvar_bg.to_csv(outdir / "chromvar_bagging_tflist.txt", sep="\t", index=True)

    chromvarBG = chromvar[:, chromvar_bg.index].copy()
    chromvarBG.write(outdir / "chromvarBG.h5ad")

    print(f"done: {chromvar.shape[1]} motifs scanned, {chromvarBG.shape[1]} retained after bagging")


if __name__ == "__main__":
    main()
