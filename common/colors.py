"""Shared helpers for deterministic default color assignment."""

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt


def default_color_dict(categories):
    """Deterministic fallback palette (tab20+tab20b+tab20c, 60 colors, cycled)
    for categories with no explicit color map."""
    colors = []
    for name in ("tab20", "tab20b", "tab20c"):
        colors.extend(mcolors.rgb2hex(c) for c in plt.get_cmap(name).colors)
    return {cat: colors[i % len(colors)] for i, cat in enumerate(sorted(categories))}
