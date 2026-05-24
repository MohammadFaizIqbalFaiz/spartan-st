# Recommended parameter settings

| Technology | `spatial_coord` | `spatial_neighborhood` | `spatial_neighs` | `spatial_rings` | `gene_coord` | `gene_neighs` | Suggested `alpha` | Suggested `beta1` | Suggested `beta2` | Notes |
|---|---|---|---:|---:|---|---:|---|---|---|---|
| Visium HD | `grid` | `knn` | 4–6 | 2 | `generic` | 15 | 0.70–0.80 | 0.10–0.26 | 0.24–0.40 | Use higher `alpha` values for large high-resolution datasets. |
| MERFISH | `generic` | `knn` | 10–12 | NA | `generic` | 15 | 0.69–0.82 | 0.10–0.26 | 0.24–0.40 | Recommended for single-cell imaging datasets. |
| MERFISH | `generic` | `delaunay` | NA | NA | `generic` | 15 | 0.75–0.85 | 0.10–0.26 | 0.24–0.40 | Useful for irregular imaging-based layouts. |
| Stereo-seq | `grid` | `knn` | 4 | 1 | `generic` | 15 | 0.50–0.60 | 0.26 | 0.24 | Sequencing-based datasets often stabilize at lower `alpha`. |
| Visium SD | `grid` | `knn` | 6 | 2 | `generic` | 15 | 0.55–0.75 | 0.26 | 0.24 | Standard Visium spot-level datasets benefit from moderate `alpha`. |

## Parameters users can tune

| Parameter | Default | Description |
|---|---:|---|
| `gene_neighs` | 15 | Number of neighbors used to construct the gene expression connectivity graph. |
| `total_pca_comps` | 50 | Total number of principal components computed during PCA. |
| `pca_comps_extract` | 30 | Number of principal components used for graph construction. |
