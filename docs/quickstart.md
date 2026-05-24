# Quickstart tutorial

This quickstart is intended for new Spartan users. It demonstrates a minimal end-to-end workflow using the representative MERFISH dataset available from Squidpy. The workflow loads one MERFISH section, preprocesses the data, identifies Spartan spatial domains, runs SAQ-based spatially variable gene discovery, and visualizes both spatial domains and representative SVGs.

---

## 1. Import packages and load a MERFISH sample

```python
import spartan as sp
import squidpy as sq
import scanpy as sc
import matplotlib.pyplot as plt

adata_data = sq.datasets.merfish()  # representative MERFISH dataset from Squidpy

# Extract one sample/section
adata = adata_data[adata_data.obs["Bregma"] == 1].copy()
```

---

## 2. Add a ground-truth column if needed

Some Spartan workflows and tutorials expect a `ground_truth` column for downstream benchmarking or visualization. The core Spartan domain workflow itself does not require ground-truth labels.

```python
# If a ground-truth column is absent, initialize one
adata.obs["ground_truth"] = ""
```

---

## 3. Preprocess the imaging-based dataset

For MERFISH and related imaging-based spatial transcriptomics datasets, use the imaging preprocessing wrapper. To use this wrapper, a `ground_truth` column is expected. If your dataset does not already contain one, create it using the code snippet above.

```python
adata = sp.tl.pre_process_imaging(adata)
```

---

## 4. Run Spartan spatial domain identification

This step constructs the spatial graph, Local Spatial Activation graph, gene-expression connectivity graph, aggregated Spartan graph, and then applies Leiden clustering to identify spatial domains.

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

The predicted spatial domains are stored in:

```python
adata.obs["spartan_domains"]
```

Inspect the number of cells assigned to each Spartan domain:

```python
adata.obs["spartan_domains"].value_counts()
```

---

## 5. Run SAQ-based spatially variable gene discovery

Spartan uses the Local Spatial Activation graph to compute Spatial Activation Quotient scores for spatially variable gene discovery.

```python
sp.tl.spartan_svg(
    adata,
    lsa_graph=adata.obsp["spartan_lsa_graph"],
    top_k_refine=161,
    key_added="spartan_svg",
)
```

View the top SAQ-ranked genes:

```python
adata.var.sort_values("spartan_saq", ascending=False).head(10)
```

---

## 6. Visualize cell classes and Spartan spatial domains

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

Example output:

```{image} _static/quickstart_domains.png
:width: 900px
:align: center
```

---

## 7. Visualize representative Spartan SVGs

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

Example output:

```{image} _static/quickstart_svgs.png
:width: 900px
:align: center
```

---

## 8. Main objects created by the quickstart

After running the quickstart, Spartan stores the main outputs inside the `AnnData` object.

### Spatial domain labels

```python
adata.obs["spartan_domains"]
```

### Graph outputs

```python
adata.obsp["spartan_spatial_graph"]
adata.obsp["spartan_spatial_weights"]
adata.obsp["spartan_lsa_graph"]
adata.obsp["spartan_gene_graph"]
adata.obsp["spartan_joint_graph"]
```

### SAQ/SVG outputs

```python
adata.var["spartan_saq"]
adata.var["spartan_saq_pval"]
adata.var["spartan_saq_fdr"]
adata.var["spartan_svg"]
adata.var["spartan_saq_rank"]
```

---

## Complete quickstart code

```python
import spartan as sp
import squidpy as sq
import scanpy as sc
import matplotlib.pyplot as plt

adata_data = sq.datasets.merfish()
adata = adata_data[adata_data.obs["Bregma"] == 1].copy()

adata.obs["ground_truth"] = ""

adata = sp.tl.pre_process_imaging(adata)

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

adata.obs["spartan_domains"].value_counts()

sp.tl.spartan_svg(
    adata,
    lsa_graph=adata.obsp["spartan_lsa_graph"],
    top_k_refine=161,
    key_added="spartan_svg",
)

adata.var.sort_values("spartan_saq", ascending=False).head(10)

sc.pl.spatial(
    adata,
    img_key=None,
    color=["Cell_class", "spartan_domains"],
    wspace=0.4,
    title=["Cell Class", "Spatial Domains"],
    size=2.5,
    spot_size=0.005,
)

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

For sequencing-based spatial transcriptomics, use:

```python
adata = sp.tl.pre_process_sequencing(adata)
```
