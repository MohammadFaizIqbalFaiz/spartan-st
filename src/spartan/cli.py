from __future__ import annotations

import argparse
from pathlib import Path
import anndata as ad

from .tl import spartan_spatial_domains, spartan_svg


def main(argv=None):
    p = argparse.ArgumentParser(prog="spartan", description="SPARTAN CLI (basic).")
    sub = p.add_subparsers(dest="cmd", required=True)

    fit = sub.add_parser("domains", help="Run spatial domains identification on an AnnData .h5ad")
    fit.add_argument("--h5ad", required=True, type=Path)
    fit.add_argument("--out", required=True, type=Path)
    fit.add_argument("--key-added", default="spartan_domains")
    fit.add_argument("--alpha", type=float, default=0.80)
    fit.add_argument("--beta1", type=float, default=0.10)
    fit.add_argument("--beta2", type=float, default=0.40)
    fit.add_argument("--resolution", type=float, default=1.0)
    fit.add_argument("--seed", type=int, default=1)

    svg = sub.add_parser("svg", help="Run SVG discovery (requires LSA graph already computed)")
    svg.add_argument("--h5ad", required=True, type=Path)
    svg.add_argument("--out", required=True, type=Path)
    svg.add_argument("--lsa-key", default="spartan_lsa_graph")
    svg.add_argument("--key-added", default="spartan_svg")

    args = p.parse_args(argv)

    if args.cmd == "domains":
        adata = ad.read_h5ad(args.h5ad)
        spartan_spatial_domains(
            adata,
            alpha=args.alpha, beta1=args.beta1, beta2=args.beta2,
            resolution=args.resolution, seed=args.seed,
            key_added=args.key_added,
            copy=False,
        )
        args.out.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(args.out)
        return 0

    if args.cmd == "svg":
        adata = ad.read_h5ad(args.h5ad)
        lsa = adata.obsp[args.lsa_key]
        spartan_svg(adata, lsa_graph=lsa, key_added=args.key_added, copy=False)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(args.out)
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
