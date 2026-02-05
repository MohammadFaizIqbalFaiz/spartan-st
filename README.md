<p align="center">
  <img src="assets/spartan-logo.svg" width="520"/>
</p>

# SPARTAN

**SPARTAN** (Spatial Activation–Aware Transcriptomic Analysis Network) is a Python toolkit
for **spatial domain identification** and **spatially variable gene (SVG) discovery**
in spatial transcriptomics data.

SPARTAN integrates spatial proximity, gene expression similarity, and **local spatial
activation (LSA)** — a neighborhood-dependent spatial autocorrelation signal —
to resolve biologically coherent tissue domains and identify genes associated with
localized spatial structure. The method is designed to preserve anatomical boundaries,
detect transition zones, and highlight spatial microenvironments across diverse spatial
transcriptomics technologies.

---

## Conceptual overview

Most spatial clustering methods rely on two signals:
(1) physical proximity and (2) transcriptomic similarity.
SPARTAN introduces a third signal — **local spatial activation (LSA)** —
which quantifies how strongly a spot’s molecular state deviates from its local neighborhood.

SPARTAN constructs three graphs:

1. **Spatial graph (S)**  
   Encodes physical proximity between neighboring spots or cells.

2. **Gene expression graph (G)**  
   Captures transcriptomic similarity in a latent (PCA) space.

3. **Local spatial activation graph (L)**  
   Measures neighborhood-dependent spatial autocorrelation, highlighting tissue
   interfaces, boundaries, and transition regions.

These graphs are combined into a **multiplex joint graph**, which is partitioned using
the Leiden algorithm to identify **spatial domains**.  
The same local activation signal is used to compute a **Spatial Activation Quotient (SAQ)**
for robust, permutation-based detection of **spatially variable genes**.

---

## Key features

- Multiplex graph-based spatial domain identification
- Explicit modeling of local spatial activation
- Robust SVG discovery using the Spatial Activation Quotient (SAQ)
- Scanpy-style API (`spartan.tl`, `spartan.pl`)
- Compatible with AnnData and SpatialData tables
- Scales across sequencing- and imaging-based spatial technologies

---

## Installation

### Recommended (conda)

```bash
conda env create -f environment.yml
conda activate biogis
pip install -e .
```
### Pip 
```bash
pip install spartan-st
```

## Quickstart (AnnData)

```python
import spartan as sp

# Spatial domain identification
sp.tl.spartan_spatial_domains(
    adata,
    key_added="spartan_domains",
)

# Inspect domain assignments
adata.obs["spartan_domains"].value_counts()

# Spatially variable gene discovery
sp.tl.spartan_svg(
    adata,
    lsa_graph=adata.obsp["spartan_lsa_graph"],
)

# View top spatially variable genes
adata.var.sort_values("spartan_saq", ascending=False).head(10)

## Outputs

### Spatial domains
- `adata.obs["spartan_domains"]` — spatial domain labels for each spot or cell

### Graphs (stored in `adata.obsp`)
- `adata.obsp["spartan_spatial_graph"]` — spatial proximity graph
- `adata.obsp["spartan_spatial_weights"]` — spatial weight matrix
- `adata.obsp["spartan_lsa_graph"]` — local spatial activation graph
- `adata.obsp["spartan_gene_graph"]` — gene expression similarity graph
- `adata.obsp["spartan_joint_graph"]` — multiplex joint graph used for clustering

### Spatially variable genes (SVGs)
- `adata.var["spartan_saq"]` — Spatial Activation Quotient (SAQ) score
- `adata.var["spartan_saq_pval"]` — permutation-based p-value
- `adata.var["spartan_saq_fdr"]` — FDR-adjusted p-value
- `adata.var["spartan_svg"]` — boolean indicator of significant SVGs
- `adata.var["spartan_saq_rank"]` — SAQ-based gene ranking

---

## Tested environment

SPARTAN was developed and tested with the following package versions:

- numpy 2.2.6
- scipy 1.15.2
- anndata 0.11.4
- scanpy 1.11.4
- squidpy 1.6.5
- igraph 0.11.8
- leidenalg 0.10.2
- joblib 1.5.1
- statsmodels 0.14.5
- matplotlib 3.10.5
- spatialdata 0.4.0

For exact reproducibility, use the provided `environment.yml`.





