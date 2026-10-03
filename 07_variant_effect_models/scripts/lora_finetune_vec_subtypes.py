#!/usr/bin/env python
"""LoRA fine-tune a per-celltype seq2PRINT base model on its subtypes.

Builds a multi-group LoRA fine-tuning config -- one embedding-conditioned
group per subtype -- on top of an already-trained per-celltype base model
(see train_seq2print_per_celltype.py), then prints (or, with --launch,
submits) the training command.

Requires a scPrinter printer object for the base cell type, its peak set,
its trained base model + config, a barcode -> subtype assignment (already
filtered to subtypes with enough sequencing depth), and a per-barcode
embedding table (e.g. RNA PCA) to condition each subtype's LoRA group --
reconciling RNA subtype annotations with ATAC barcodes and filtering by
depth is a separate, cohort-specific data-prep step and out of scope here.

Consolidated from seq2print/1-1b2_base_model_perCelltype.ipynb sections
3-4 (the LoRA fine-tune step, 4-2-1). Originally used to fine-tune the
Blood Vascular Endothelial Cell base model on its VEC subtypes, but the
method is general to any base cell type / subtype grouping.
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scprinter as scp


def sanitize_subtype(subtype: str) -> str:
    return re.sub(r"[^\w]", "_", subtype.replace(" ", "_"))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--printer", required=True, help="scPrinter printer h5ad for the base cell type")
    ap.add_argument("--region-path", required=True, help="peak set BED used to train the base model (same cell type)")
    ap.add_argument("--pretrain-model", required=True, help="trained base model .pt (see train_seq2print_per_celltype.py)")
    ap.add_argument("--base-model-config", required=True, help="the base model's training config JSON")
    ap.add_argument(
        "--subtype-barcodes",
        required=True,
        help="tab-separated table with columns 'barcode', 'subtype' (already depth-filtered to subtypes worth fine-tuning on)",
    )
    ap.add_argument(
        "--embeddings",
        required=True,
        help="tab-separated table, barcode in the first column, embedding dims (e.g. RNA PCA) as the remaining columns",
    )
    ap.add_argument("--genome", default="hg38", choices=["hg38", "mm10"], help="scprinter.genome build (default: hg38)")
    ap.add_argument("--work-dir", required=True, help="seq2PRINT project directory; configs/model/temp are created under it")
    ap.add_argument("--model-name", default="CZI_LoRA", help="LoRA model name, used in the config filename and wandb tags (default: CZI_LoRA)")
    ap.add_argument("--n-folds", type=int, default=1, help="number of training folds (default: 1)")
    ap.add_argument("--lr", type=float, default=3e-5, help="LoRA learning rate (default: 3e-5)")
    ap.add_argument("--gpus", type=int, default=1, help="number of GPUs for the training job (default: 1)")
    ap.add_argument("--wandb-project", help="wandb project name (default: scPrinter_{model-name})")
    ap.add_argument(
        "--launch",
        action="store_true",
        help="actually submit the training job instead of just printing the command (default: off -- this is an "
        "expensive, long-running GPU job)",
    )
    args = ap.parse_args()

    work_dir = Path(args.work_dir)
    configs_dir = work_dir / "configs"
    configs_dir.mkdir(parents=True, exist_ok=True)

    genome = getattr(scp.genome, args.genome)
    printer = scp.load_printer(args.printer, genome)

    subtype_barcodes = pd.read_csv(args.subtype_barcodes, sep="\t")
    subtype_barcodes["subtype_clean"] = subtype_barcodes["subtype"].map(sanitize_subtype)

    embeddings_all = pd.read_csv(args.embeddings, sep="\t", index_col=0)

    cell_grouping_all, group_names_all = [], []
    for subtype_org, subtype_clean in (
        subtype_barcodes[["subtype", "subtype_clean"]].drop_duplicates().itertuples(index=False)
    ):
        barcodes = subtype_barcodes.loc[subtype_barcodes["subtype"] == subtype_org, "barcode"].values
        barcodes_in_embeddings = np.intersect1d(barcodes, embeddings_all.index.values)
        print(f"{subtype_org} (clean: {subtype_clean}): {len(barcodes)} barcodes, {len(barcodes_in_embeddings)} with embeddings")
        cell_grouping_all.append(barcodes_in_embeddings)
        group_names_all.append(subtype_clean)

    all_barcodes_used = np.concatenate(cell_grouping_all)
    embeddings_sub = embeddings_all.loc[np.unique(all_barcodes_used)]

    wandb_project = args.wandb_project or f"scPrinter_{args.model_name}"

    for fold in range(args.n_folds):
        config_path = configs_dir / f"{args.model_name}.fold{fold}.JSON"
        scp.tl.seq_lora_model_config(
            printer,
            region_path=args.region_path,
            cell_grouping=cell_grouping_all,
            group_names=group_names_all,
            embeddings=embeddings_sub,
            genome=genome,
            pretrain_model=args.pretrain_model,
            overwrite_barcode=False,
            model_name=args.model_name,
            fold=fold,
            model_config=args.base_model_config,
            additional_lora_config={
                "lr": args.lr,
                "notes": "multi-group LoRA model",
                "tags": ["CZI", "perCelltype", "LoRA", "multi-group", args.model_name, f"fold{fold}"],
            },
            config_save_path=str(config_path),
        )

        scp.tl.launch_seq2print(
            model_config_path=str(config_path),
            temp_dir=str(work_dir / "temp"),
            model_dir=str(work_dir / "model"),
            data_dir=str(work_dir),
            gpus=args.gpus,
            wandb_project=wandb_project,
            verbose=False,
            launch=args.launch,
        )

    printer.close()
    print("done")


if __name__ == "__main__":
    main()
