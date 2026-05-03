"""Command line interface for SPARTAN."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad

from .tl import (
    pre_process_imaging,
    pre_process_sequencing,
    spartan_build_graphs,
    spartan_spatial_domains,
    spartan_svg,
    initiate_alpha_selection,
)


def _add_graph_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--spatial-coord", default="grid")
    parser.add_argument("--spatial-neighs", type=int, default=6)
    parser.add_argument("--spatial-rings", type=int, default=2)
    parser.add_argument("--spatial-neighborhood", default="knn", choices=["knn", "delaunay"])
    parser.add_argument("--total-pca-comps", type=int, default=50)
    parser.add_argument("--pca-comps-extract", type=int, default=30)
    parser.add_argument("--gene-coord", default="generic")
    parser.add_argument("--gene-neighs", type=int, default=15)
    parser.add_argument("--seed", type=int, default=1)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="spartan", description="SPARTAN CLI.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    pp_img = sub.add_parser("preprocess-imaging", help="Preprocess an imaging-based SRT AnnData .h5ad.")
    pp_img.add_argument("--h5ad", required=True, type=Path)
    pp_img.add_argument("--out", required=True, type=Path)

    pp_seq = sub.add_parser("preprocess-sequencing", help="Preprocess a sequencing-based SRT AnnData .h5ad.")
    pp_seq.add_argument("--h5ad", required=True, type=Path)
    pp_seq.add_argument("--out", required=True, type=Path)
    pp_seq.add_argument("--min-counts", type=int, default=500)
    pp_seq.add_argument("--min-cells", type=int, default=5)
    pp_seq.add_argument("--jitter", type=float, default=0.4)
    pp_seq.add_argument("--n-top-genes", type=int, default=3000)
    pp_seq.add_argument("--seed", type=int, default=1)

    build = sub.add_parser("build-graphs", help="Build SPARTAN spatial, LSA, and gene graphs.")
    build.add_argument("--h5ad", required=True, type=Path)
    build.add_argument("--out", required=True, type=Path)
    _add_graph_args(build)

    domains = sub.add_parser("domains", help="Run SPARTAN spatial domain identification.")
    domains.add_argument("--h5ad", required=True, type=Path)
    domains.add_argument("--out", required=True, type=Path)
    domains.add_argument("--key-added", default="spartan_domains")
    domains.add_argument("--alpha", type=float, default=0.80)
    domains.add_argument("--beta1", type=float, default=0.10)
    domains.add_argument("--beta2", type=float, default=0.40)
    domains.add_argument("--resolution", type=float, default=1.0)
    _add_graph_args(domains)

    svg = sub.add_parser("svg", help="Run SPARTAN SAQ/SVG discovery.")
    svg.add_argument("--h5ad", required=True, type=Path)
    svg.add_argument("--out", required=True, type=Path)
    svg.add_argument("--lsa-key", default="spartan_lsa_graph")
    svg.add_argument("--layer", default="log1pX")
    svg.add_argument("--key-added", default="spartan_svg")
    svg.add_argument("--n-permutations", type=int, default=1000)
    svg.add_argument("--n-cores", type=int, default=8)
    svg.add_argument("--alpha-svg", type=float, default=0.05)
    svg.add_argument("--chunk-size", type=int, default=200)
    svg.add_argument("--seed", type=int, default=1)
    svg.add_argument("--single-stage", action="store_true", help="Disable two-stage SVG refinement.")
    svg.add_argument("--n-permutations-stage1", type=int, default=100)
    svg.add_argument("--top-k-refine", type=int, default=3000)
    svg.add_argument("--prefer-backend", default="threads", choices=["threads", "processes"])

    alpha = sub.add_parser("alpha-select", help="Run SPARTAN alpha-selection on a graph-built AnnData.")
    alpha.add_argument("--h5ad", required=True, type=Path)
    alpha.add_argument("--out-csv", required=True, type=Path)
    alpha.add_argument("--lower-alpha", type=float, required=True)
    alpha.add_argument("--upper-alpha", type=float, required=True)
    alpha.add_argument("--step-alpha", type=float, required=True)
    alpha.add_argument("--lower-resolution", type=float, required=True)
    alpha.add_argument("--upper-resolution", type=float, required=True)
    alpha.add_argument("--step-resolution", type=float, required=True)
    alpha.add_argument("--lower-nlsas", type=int, default=40)
    alpha.add_argument("--upper-nlsas", type=int, default=70)
    alpha.add_argument("--n-jobs", type=int, default=4)
    alpha.add_argument("--config", default="lsg", choices=["lsg", "lg", "sg"])
    alpha.add_argument("--seed", type=int, default=1)
    alpha.add_argument("--use-nLSAS", action="store_true")

    args = parser.parse_args(argv)

    if args.cmd == "preprocess-imaging":
        adata = ad.read_h5ad(args.h5ad)
        adata = pre_process_imaging(adata)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(args.out)
        return 0

    if args.cmd == "preprocess-sequencing":
        adata = ad.read_h5ad(args.h5ad)
        adata = pre_process_sequencing(
            adata,
            min_counts=args.min_counts,
            min_cells=args.min_cells,
            jitter=args.jitter,
            n_top_genes=args.n_top_genes,
            seed=args.seed,
        )
        args.out.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(args.out)
        return 0

    if args.cmd == "build-graphs":
        adata = ad.read_h5ad(args.h5ad)
        spartan_build_graphs(
            adata,
            spatial_coord=args.spatial_coord,
            spatial_neighs=args.spatial_neighs,
            spatial_rings=args.spatial_rings,
            spatial_neighborhood=args.spatial_neighborhood,
            total_pca_comps=args.total_pca_comps,
            pca_comps_extract=args.pca_comps_extract,
            gene_coord=args.gene_coord,
            gene_neighs=args.gene_neighs,
            seed=args.seed,
            copy=False,
        )
        args.out.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(args.out)
        return 0

    if args.cmd == "domains":
        adata = ad.read_h5ad(args.h5ad)
        spartan_spatial_domains(
            adata,
            spatial_coord=args.spatial_coord,
            spatial_neighs=args.spatial_neighs,
            spatial_rings=args.spatial_rings,
            spatial_neighborhood=args.spatial_neighborhood,
            total_pca_comps=args.total_pca_comps,
            pca_comps_extract=args.pca_comps_extract,
            gene_coord=args.gene_coord,
            gene_neighs=args.gene_neighs,
            alpha=args.alpha,
            beta1=args.beta1,
            beta2=args.beta2,
            resolution=args.resolution,
            seed=args.seed,
            key_added=args.key_added,
            copy=False,
        )
        args.out.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(args.out)
        return 0

    if args.cmd == "svg":
        adata = ad.read_h5ad(args.h5ad)
        spartan_svg(
            adata,
            lsa_graph=adata.obsp[args.lsa_key],
            layer=args.layer,
            n_permutations=args.n_permutations,
            n_cores=args.n_cores,
            alpha_svg=args.alpha_svg,
            chunk_size=args.chunk_size,
            seed=args.seed,
            key_added=args.key_added,
            two_stage=not args.single_stage,
            n_permutations_stage1=args.n_permutations_stage1,
            top_k_refine=args.top_k_refine,
            prefer_backend=args.prefer_backend,
            copy=False,
        )
        args.out.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(args.out)
        return 0

    if args.cmd == "alpha-select":
        adata = ad.read_h5ad(args.h5ad)
        summary = initiate_alpha_selection(
            adata,
            lower_alpha=args.lower_alpha,
            upper_alpha=args.upper_alpha,
            lower_resolution=args.lower_resolution,
            upper_resolution=args.upper_resolution,
            step_alpha=args.step_alpha,
            step_resolution=args.step_resolution,
            lower_nlsas=args.lower_nlsas,
            upper_nlsas=args.upper_nlsas,
            n_jobs=args.n_jobs,
            config=args.config,
            seed=args.seed,
            use_nLSAS=args.use_nLSAS,
        )
        args.out_csv.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(args.out_csv, index=False)
        return 0

    raise RuntimeError(f"Unknown command: {args.cmd}")


if __name__ == "__main__":
    raise SystemExit(main())
