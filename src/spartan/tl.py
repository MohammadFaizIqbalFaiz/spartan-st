"""Tool functions (Scanpy-like)."""

from __future__ import annotations
from typing import Optional
from anndata import AnnData
from scipy.sparse import csr_matrix

from ._core import (
    spartan_spatial_domains,
    spartan_svg,
)

__all__ = ["spartan_spatial_domains", "spartan_svg"]
