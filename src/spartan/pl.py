"""Plotting helpers for SPARTAN.

The plotting API is intentionally lightweight and Scanpy-compatible.
"""

from __future__ import annotations

from typing import Optional
import matplotlib.pyplot as plt
from anndata import AnnData


def spatial_domains(
    adata: AnnData,
    color: str = "spartan_domains",
    ax=None,
    title: Optional[str] = None,
    **kwargs,
):
    """Plot Spartan spatial domains using :func:`scanpy.pl.spatial`.

    Requires ``adata.obsm["spatial"]`` to contain spatial coordinates.
    """
    import scanpy as sc

    if ax is None:
        ax = plt.gca()

    sc.pl.spatial(adata, color=color, ax=ax, title=title, show=False, **kwargs)
    return ax


def local_spatial_activation(
    adata: AnnData,
    color: str = "spartan_local_spatial_activation",
    ax=None,
    title: Optional[str] = "Local Spatial Activation",
    cmap: str = "Reds",
    **kwargs,
):
    """Plot per-spot/cell Local Spatial Activation values."""
    import scanpy as sc

    if ax is None:
        ax = plt.gca()

    sc.pl.spatial(
        adata,
        color=color,
        ax=ax,
        title=title,
        cmap=cmap,
        show=False,
        **kwargs,
    )
    return ax


def svg_table(adata: AnnData, n: int = 20, score_key: str = "spartan_saq"):
    """Return top-n SVGs as a DataFrame sorted by SAQ score."""
    if score_key not in adata.var:
        raise KeyError(f"{score_key!r} not found in adata.var.")
    return adata.var.sort_values(score_key, ascending=False).head(n)


def alpha_selection_summary(summary_df, metric: str = "median_NMI", ax=None):
    """Plot alpha-selection summary from ``spartan.tl.initiate_alpha_selection``."""
    if metric not in summary_df:
        raise KeyError(f"{metric!r} not found in summary dataframe.")

    if ax is None:
        ax = plt.gca()

    ax.plot(summary_df["alpha"], summary_df[metric], marker="o")
    ax.set_xlabel("alpha")
    ax.set_ylabel(metric)
    ax.set_title("SPARTAN alpha-selection summary")
    return ax
