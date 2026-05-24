# Spatially variable gene discovery

Spartan identifies spatially variable genes using the Spatial Activation Quotient
(SAQ). SAQ measures how strongly each gene expression pattern aligns with the
Local Spatial Activation graph.

```python
sp.tl.spartan_svg(
    adata,
    lsa_graph=adata.obsp["spartan_lsa_graph"],
    layer="log1pX",
    n_permutations=1000,
    n_cores=8,
    use_X_if_missing=True,
    alpha_svg=0.05,
    chunk_size=200,
    seed=1,
    key_added="spartan_svg",
    copy=False,
    dtype=np.float32,
    prefer_backend="threads",
    two_stage=True,
    n_permutations_stage1=100,
    top_k_refine=3000,
)
```

---

## SAQ score

For each gene \(j\), Spartan computes:

```{math}
SAQ_j = \frac{x_j'^T L x_j'}{\lVert x_j' \rVert^2}
```

where:

- \(x_j'\) is the mean-centered expression vector of gene \(j\),
- \(L\) is the Local Spatial Activation graph.

Higher SAQ values indicate stronger alignment between the gene expression pattern
and the activation-aware spatial graph.

---

## Key parameters

| Parameter | Default | Description |
|---|---:|---|
| `adata` | required | Input `AnnData` object containing expression values and gene metadata. |
| `lsa_graph` | required | Local Spatial Activation graph, usually `adata.obsp["spartan_lsa_graph"]`, used to score spatial activation of each gene. |
| `layer` | `"log1pX"` | Expression layer used for SAQ scoring. If this layer is unavailable and `use_X_if_missing=True`, Spartan falls back to `adata.X`. |
| `n_permutations` | `1000` | Total number of permutations used to estimate the null distribution of SAQ scores. Higher values provide more stable p-values but increase runtime. |
| `n_cores` | `8` | Number of CPU cores used for parallel permutation testing. |
| `use_X_if_missing` | `True` | If `True`, uses `adata.X` when the specified layer is not found. If `False`, raises an error when the layer is missing. |
| `alpha_svg` | `0.05` | FDR threshold used to call significant SVGs. Genes with FDR below this threshold are marked as spatially variable. |
| `chunk_size` | `200` | Number of genes processed per chunk. Smaller chunks reduce memory use; larger chunks may improve speed on machines with more memory. |
| `seed` | `1` | Random seed used for reproducible permutation testing. |
| `key_added` | `"spartan_svg"` | Column name added to `adata.var` containing Boolean SVG calls. |
| `copy` | `False` | If `True`, returns a modified copy of `adata`. If `False`, updates `adata` in place. |
| `dtype` | `np.float32` | Dense array type used during chunked computation. `float32` reduces memory usage and is usually sufficient. |
| `prefer_backend` | `"threads"` | Joblib backend used for parallel permutation testing. Use `"threads"` for lower overhead or `"processes"` if thread contention occurs. |
| `two_stage` | `True` | Enables two-stage permutation testing. A small first-stage run is applied to all genes, followed by refinement of top candidates. |
| `n_permutations_stage1` | `100` | Number of first-stage permutations applied to all genes when `two_stage=True`. |
| `top_k_refine` | `3000` | Number of top candidate genes refined using the remaining permutations in stage 2. |

---

## Outputs

After running `sp.tl.spartan_svg`, Spartan stores the following columns in
`adata.var`:

| Output column | Description |
|---|---|
| `spartan_saq` | Observed SAQ score for each gene. |
| `spartan_saq_pval` | One-sided permutation-based p-value for each gene. |
| `spartan_saq_fdr` | Benjamini-Hochberg adjusted FDR value. |
| `spartan_svg` or `key_added` | Boolean SVG call using the selected `alpha_svg` threshold. |
| `spartan_saq_rank` | Rank of each gene by decreasing SAQ score. |

---

## Recommended usage

```python
sp.tl.spartan_svg(
    adata,
    lsa_graph=adata.obsp["spartan_lsa_graph"],
    layer="log1pX",
    n_permutations=1000,
    n_cores=8,
    alpha_svg=0.05,
    key_added="spartan_svg",
)
```

For a faster exploratory run:

```python
sp.tl.spartan_svg(
    adata,
    lsa_graph=adata.obsp["spartan_lsa_graph"],
    n_permutations=200,
    n_permutations_stage1=50,
    top_k_refine=500,
    n_cores=8,
)
```

For a more stringent final analysis:

```python
sp.tl.spartan_svg(
    adata,
    lsa_graph=adata.obsp["spartan_lsa_graph"],
    n_permutations=5000,
    n_permutations_stage1=500,
    top_k_refine=3000,
    n_cores=16,
)
```

Ground-truth labels are not used during SAQ/SVG discovery.
