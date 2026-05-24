# Spatially variable gene discovery

Spartan identifies spatially variable genes using the Spatial Activation Quotient (SAQ).

```python
sp.tl.spartan_svg(
    adata,
    lsa_graph=adata.obsp["spartan_lsa_graph"],
    layer="log1pX",
    n_permutations=1000,
    n_cores=8,
    alpha_svg=0.05,
    chunk_size=200,
    seed=1,
    key_added="spartan_svg",
    copy=False,
)
```

## Key parameters

| Parameter | Default | Description |
|---|---:|---|
| `adata` | required | Input `AnnData` object. |
| `lsa_graph` | required | Local Spatial Activation graph used to score spatial activation of each gene. |
| `layer` | `"log1pX"` | Expression layer used for SAQ scoring. |
| `n_permutations` | `1000` | Total permutations used to estimate the null distribution. |
| `n_cores` | `8` | CPU cores used for parallel permutation testing. |
| `alpha_svg` | `0.05` | FDR threshold used to call significant SVGs. |
| `chunk_size` | `200` | Number of genes processed per chunk. |
| `two_stage` | `True` | Enables two-stage permutation testing. |
| `n_permutations_stage1` | `100` | First-stage permutations applied to all genes. |
| `top_k_refine` | `3000` | Top candidate genes refined in stage 2. |

## Outputs

- `adata.var["spartan_saq"]`
- `adata.var["spartan_saq_pval"]`
- `adata.var["spartan_saq_fdr"]`
- `adata.var["spartan_svg"]` or `adata.var[key_added]`
- `adata.var["spartan_saq_rank"]`
