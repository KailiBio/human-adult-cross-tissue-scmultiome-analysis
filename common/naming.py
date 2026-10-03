"""Shared helpers for naming output files consistently across scripts."""

import re


def sanitize_label(label: str) -> str:
    """Make a free-text label (e.g. a celltype or subtype name) filesystem-safe.

    Collapses any run of non-alphanumeric/underscore characters to a single
    "_" and strips leading/trailing "_", e.g. "Blood Vascular Endothelial
    Cell" -> "Blood_Vascular_Endothelial_Cell", "T/NK cell" -> "T_NK_cell".
    """
    return re.sub(r"[^\w]+", "_", label).strip("_")
