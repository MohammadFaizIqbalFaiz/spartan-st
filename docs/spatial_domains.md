# Spatial domain identification

The main user-facing function for spatial domain identification is:

```python
sp.tl.spartan_spatial_domains(
    adata,
    spatial_coord="grid",
    spatial_neighs=6,
    spatial_rings=2,
    spatial_neighborhood="knn",
    total_pca_comps=50,
    pca_comps_extract=30,
    gene_coord="generic",
    gene_neighs=15,
    alpha=0.80,
    beta1=0.10,
    beta2=0.40,
    resolution=1.0,
    seed=1,
    key_added="spartan_domains",
    copy=False,
)
```

`spartan_spatial_domains` is the high-level wrapper that constructs the Spartan
graph layers and performs Leiden clustering. It updates the input `AnnData`
object in place by default.

---

## Workflow

The function performs the following steps:

1. Builds the spatial graph `S`.
2. Builds the row-normalized spatial weight matrix.
3. Performs PCA on the preprocessed expression matrix.
4. Builds the Local Spatial Activation graph `L`.
5. Builds the gene expression connectivity graph `G`.
6. Aggregates the graphs into the Spartan joint graph `J`.
7. Runs Leiden clustering on `J`.
8. Stores the resulting spatial domains in `adata.obs[key_added]`.

The aggregated graph is:

```{math}
J = (\alpha-\beta_1)L + (1-\alpha)G + (\alpha-\beta_2)S
```

where:

- `L` is the Local Spatial Activation graph,
- `G` is the gene expression connectivity graph,
- `S` is the spatial adjacency graph.

The effective graph contributions are therefore:

```python
LSA contribution     = alpha - beta1
Gene contribution    = 1 - alpha
Spatial contribution = alpha - beta2
```

---

## Key parameters

| Parameter | Default | Description |
|---|---:|---|
| `adata` | required | Input `AnnData` object. It should contain expression values and spatial coordinates in `adata.obsm["spatial"]`. |
| `spatial_coord` | `"grid"` | Coordinate type used for spatial graph construction. Use `"grid"` for Visium, Visium HD, and Stereo-seq-style gridded data; use `"generic"` for imaging-based single-cell data such as MERFISH or Vizgen MERFISH. |
| `spatial_neighs` | `6` | Number of spatial neighbors used when `spatial_neighborhood="knn"`. Larger values produce smoother spatial connectivity; smaller values preserve more local structure. |
| `spatial_rings` | `2` | Number of spatial rings used for grid-like datasets. This is mainly relevant when `spatial_coord="grid"`. |
| `spatial_neighborhood` | `"knn"` | Spatial graph construction method. Supported values are `"knn"` and `"delaunay"`. KNN gives controlled neighborhood degree; Delaunay can be useful for irregular single-cell imaging layouts. |
| `total_pca_comps` | `50` | Total number of principal components computed from the preprocessed expression matrix. |
| `pca_comps_extract` | `30` | Number of principal components used for Spartan graph construction. |
| `gene_coord` | `"generic"` | Coordinate mode used internally when constructing the gene expression connectivity graph in PCA space. In most cases, keep this as `"generic"`. |
| `gene_neighs` | `15` | Number of neighbors used to construct the gene expression connectivity graph `G`. |
| `alpha` | `0.80` | Main graph-integration parameter controlling the balance between activation/spatial structure and gene expression connectivity. |
| `beta1` | `0.10` | Offset controlling the effective contribution of the LSA graph: `alpha - beta1`. Lower `beta1` increases the effective LSA contribution. |
| `beta2` | `0.40` | Offset controlling the effective contribution of the spatial graph: `alpha - beta2`. |
| `resolution` | `1.0` | Leiden resolution parameter. Higher values generally produce more spatial domains; lower values produce coarser domains. |
| `seed` | `1` | Random seed used for reproducibility. |
| `key_added` | `"spartan_domains"` | Column name in `adata.obs` where Spartan domain labels are stored. |
| `copy` | `False` | If `True`, returns a modified copy of `adata`. If `False`, updates `adata` in place and returns `None`. |

---

## Outputs

After running `sp.tl.spartan_spatial_domains`, Spartan stores predicted domains
and graph outputs in the `AnnData` object.

### Spatial domain labels

```python
adata.obs["spartan_domains"]
```

or:

```python
adata.obs[key_added]
```

### Graph outputs

```python
adata.obsp["spartan_spatial_graph"]
adata.obsp["spartan_spatial_weights"]
adata.obsp["spartan_lsa_graph"]
adata.obsp["spartan_gene_graph"]
adata.obsp["spartan_joint_graph"]
```

| Output | Description |
|---|---|
| `spartan_spatial_graph` | Spatial neighborhood adjacency graph `S`. |
| `spartan_spatial_weights` | Row-normalized spatial weight matrix used for LSA construction. |
| `spartan_lsa_graph` | Local Spatial Activation graph `L`. |
| `spartan_gene_graph` | Gene expression connectivity graph `G`. |
| `spartan_joint_graph` | Aggregated graph `J` used for Leiden clustering. |

---

## Recommended starting settings

### Imaging-based single-cell data, such as MERFISH or Vizgen MERFISH

```python
sp.tl.spartan_spatial_domains(
    adata,
    spatial_coord="generic",
    spatial_neighborhood="knn",
    spatial_neighs=10,
    gene_neighs=15,
    alpha=0.80,
    beta1=0.10,
    beta2=0.40,
    resolution=1.0,
    seed=1,
)
```

### Visium, Visium HD, or Stereo-seq-style grid data

```python
sp.tl.spartan_spatial_domains(
    adata,
    spatial_coord="grid",
    spatial_neighborhood="knn",
    spatial_neighs=6,
    spatial_rings=2,
    gene_neighs=15,
    alpha=0.70,
    beta1=0.26,
    beta2=0.24,
    resolution=1.0,
    seed=1,
)
```

---

## Notes

- Spartan spatial domain detection is unsupervised.
- Ground-truth labels are not used during graph construction, graph integration, or Leiden clustering.
- `alpha`, `beta1`, and `beta2` control graph integration, while `resolution` controls clustering granularity.
- Users can adjust `gene_neighs`, `total_pca_comps`, and `pca_comps_extract` depending on dataset size and biological complexity.
