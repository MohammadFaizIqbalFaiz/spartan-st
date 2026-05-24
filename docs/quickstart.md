# Quickstart tutorial

This quickstart demonstrates a minimal end-to-end workflow using the representative MERFISH dataset available from Squidpy.

## 1. Import packages and load data

```python
import spartan as sp
import squidpy as sq
import scanpy as sc
import matplotlib.pyplot as plt

adata_data = sq.datasets.merfish()
adata = adata_data[adata_data.obs["Bregma"] == 1].copy()
```

## 2. Add a ground-truth column if needed

```python
adata.obs["ground_truth"] = ""
```

## 3. Preprocess the imaging-based dataset

```python
adata = sp.tl.pre_process_imaging(adata)
```

## 4. Run Spartan spatial domain identification

```python
sp.tl.spartan_spatial_domains(
    adata,
    spatial_coord="generic",
    spatial_neighborhood="knn",
    spatial_neighs=10,
    gene_neighs=15,
    alpha=0.69,
    beta1=0.26,
    beta2=0.24,
    resolution=0.61,
    seed=1,
    key_added="spartan_domains",
)
```

## 5. Run SAQ-based SVG discovery

```python
sp.tl.spartan_svg(
    adata,
    lsa_graph=adata.obsp["spartan_lsa_graph"],
    top_k_refine=161,
    key_added="spartan_svg",
)
```

## 6. Visualize spatial domains

```python
sc.pl.spatial(
    adata,
    img_key=None,
    color=["Cell_class", "spartan_domains"],
    wspace=0.4,
    title=["Cell Class", "Spatial Domains"],
    size=2.5,
    spot_size=0.005,
)
```

```{image} _static/quickstart_domains.png
:width: 900px
:align: center
```

## 7. Visualize representative SVGs

```python
genes = ["Ucn3", "Mbp", "Sln", "Nnat"]

sc.pl.spatial(
    adata,
    img_key=None,
    color=genes,
    color_map="Reds",
    size=2.5,
    spot_size=0.005,
    legend_loc=None,
)
```

```{image} _static/quickstart_svgs.png
:width: 900px
:align: center
```
