# Spatial domain identification

The main entry point is:

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

This workflow constructs `S`, `L`, and `G`, forms the aggregated graph `J`, runs Leiden clustering, and stores spatial domains in `adata.obs[key_added]`.

Graph outputs are stored in `adata.obsp`:

- `spartan_spatial_graph`
- `spartan_spatial_weights`
- `spartan_lsa_graph`
- `spartan_gene_graph`
- `spartan_joint_graph`
