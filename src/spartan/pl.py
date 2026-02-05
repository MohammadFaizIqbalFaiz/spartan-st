"""Plotting helpers (lightweight wrappers).

For SpatialData rendering, users can use `spatialdata-plot` directly.
We keep Spartan plotting minimal and Matplotlib-based to avoid heavy deps.
"""

from __future__ import annotations
from typing import Sequence, Optional
import matplotlib.pyplot as plt
from anndata import AnnData


def spatial_domains(adata: AnnData, color: str = "spartan_domains", ax=None, title: Optional[str]=None):
    """Quick visual inspection of domains via Scanpy's spatial plot.

    Requires `adata.obsm['spatial']` to contain spatial coordinates.
    """
    import scanpy as sc
    if ax is None:
        ax = plt.gca()
    sc.pl.spatial(adata, color=color, ax=ax, title=title, show=False)
    return ax


def svg_table(adata: AnnData, n: int = 20, score_key: str = "spartan_saq"):
    """Return top-n SVGs as a DataFrame (sorted)."""
    return adata.var.sort_values(score_key, ascending=False).head(n)
