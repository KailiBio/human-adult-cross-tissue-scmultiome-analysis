#!/usr/bin/env python
"""Generate seq2PRINT training configs and launch commands, one base model
per cell type.

For each cell type in `--barcode-csv`'s `--celltype-key` column, generates
a seq2PRINT model training config (`scp.tl.seq_model_config`) using all of
that cell type's cells as one group and its own per-celltype peak set,
then prints (or, with `--launch`, submits) the training command for each
config (`scp.tl.launch_seq2print`).

Requires a scPrinter printer object and a seq2PRINT-preset peak set
already built per cell type (see sections 0-1 of
1-1b2_base_model_perCelltype.ipynb: per-celltype fragment extraction,
`scp.pp.import_fragments`, `scp.pp.call_peaks`) -- a separate, heavier
data-prep step and out of scope here.

Consolidated from seq2print/1-1b2_base_model_perCelltype.ipynb section 2.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scprinter as scp

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common.naming import sanitize_label


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--barcode-csv", required=True, help="tab-separated table with a cell type column (index=barcode)")
    ap.add_argument("--celltype-key", default="celltype", help="column in --barcode-csv with cell type labels (default: celltype)")
    ap.add_argument(
        "--work-dir",
        required=True,
        help="seq2PRINT project directory: expects 'celltype_printer/' and 'celltype_peak/' already populated per "
        "cell type; 'configs/'/'model/'/'temp/' are created under it",
    )
    ap.add_argument("--genome", default="hg38", choices=["hg38", "mm10"], help="scprinter.genome build (default: hg38)")
    ap.add_argument("--n-folds", type=int, default=1, help="number of training folds per cell type (default: 1)")
    ap.add_argument("--gpus", type=int, default=1, help="number of GPUs per training job (default: 1)")
    ap.add_argument(
        "--wandb-project-prefix",
        default="scPrinter_CZI_base",
        help="wandb project name prefix; the cell type is appended (default: scPrinter_CZI_base)",
    )
    ap.add_argument(
        "--launch",
        action="store_true",
        help="actually submit each training job instead of just printing the command (default: off -- these are "
        "expensive, long-running GPU jobs)",
    )
    ap.add_argument("--celltypes", help="comma-separated subset of cell types to process (default: all found in --barcode-csv)")
    args = ap.parse_args()

    work_dir = Path(args.work_dir)
    configs_dir = work_dir / "configs"
    configs_dir.mkdir(parents=True, exist_ok=True)

    genome = getattr(scp.genome, args.genome)

    df_barcode = pd.read_csv(args.barcode_csv, sep="\t", index_col=0, low_memory=False)
    celltypes = df_barcode[args.celltype_key].unique().tolist()
    if args.celltypes:
        wanted = set(args.celltypes.split(","))
        celltypes = [c for c in celltypes if c in wanted]
    print(f"generating seq2PRINT configs for {len(celltypes)} cell types")

    for celltype in celltypes:
        celltype_std = sanitize_label(celltype)
        print(f"processing {celltype} ({celltype_std})")

        printer_path = work_dir / "celltype_printer" / f"CZI_snATAC_base_scprinter.{celltype_std}.h5ad"
        region_path = work_dir / "celltype_peak" / f"seq2print_cleaned_narrowPeak.{celltype_std}.bed"
        if not printer_path.exists() or not region_path.exists():
            print(f"  skipping {celltype}: missing printer or peak set ({printer_path}, {region_path})")
            continue

        printer = scp.load_printer(str(printer_path), genome)
        barcode_list = np.asarray(printer.obs_names)

        for fold in range(args.n_folds):
            config_path = configs_dir / f"CZI_base.{celltype_std}.fold{fold}.JSON"
            scp.tl.seq_model_config(
                printer,
                region_path=str(region_path),
                cell_grouping=barcode_list,
                group_names=f"CZI_base.{celltype_std}",
                genome=printer.genome,
                fold=fold,
                overwrite_bigwig=False,
                model_name=f"CZI_base.{celltype_std}",
                additional_config={"notes": "v2", "tags": ["CZI", "perCelltype", celltype_std, f"fold{fold}"]},
                config_save_path=str(config_path),
            )

            scp.tl.launch_seq2print(
                model_config_path=str(config_path),
                temp_dir=str(work_dir / "temp"),
                model_dir=str(work_dir / "model"),
                data_dir=str(work_dir),
                gpus=args.gpus,
                wandb_project=f"{args.wandb_project_prefix}.{celltype_std}",
                verbose=False,
                launch=args.launch,
            )

        printer.close()

    print("done")


if __name__ == "__main__":
    main()
