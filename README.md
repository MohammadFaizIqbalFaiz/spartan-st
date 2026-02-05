# SPARTAN

**SPARTAN** provides:
- **Spatial domains identification** (graph-based clustering using joint LSA / gene / spatial graphs)
- **Spatially variable gene discovery** using the **Spatial Activation Quotient (SAQ)**

This repository is structured in a Scanpy-like way:
- `spartan.tl.*` for tools (compute + write results into `adata`)
- `spartan.pl.*` for plotting

> Note on naming: the PyPI name `spartan` is already taken by other projects, so the install name here is **`spartan-st`**. You still import it as `import spartan`.

## Install

### Conda (recommended)
```bash
conda env create -f environment.yml
conda activate biogis
pip install -e .
```

### Pip
```bash
pip install -e .
```

## Quickstart (AnnData)
```python
import spartan as sp

# domains
sp.tl.spartan_spatial_domains(adata, key_added="spartan_domains")

# inspect
adata.obs["spartan_domains"].value_counts()
adata.obsp["spartan_joint_graph"]

# SVGs
lsa = adata.obsp["spartan_lsa_graph"]
sp.tl.spartan_svg(adata, lsa_graph=lsa, key_added="spartan_svg")

# top genes
top = adata.var.sort_values("spartan_saq", ascending=False).head(20)
print(top.index.tolist())
```

## Outputs

### Domains
- `adata.obs[key_added]` : domain labels (categorical)
- Graphs in `adata.obsp`:
  - `spartan_spatial_graph`, `spartan_spatial_weights`
  - `spartan_lsa_graph`, `spartan_gene_graph`, `spartan_joint_graph`

### SVGs
Per gene in `adata.var`:
- `spartan_saq`, `spartan_saq_pval`, `spartan_saq_fdr`
- `spartan_saq_rank` (1 = highest SAQ)
- `adata.var[key_added]` boolean mask (significant SVGs)

## Citation
Add a `CITATION.cff` (template included) once you have a preprint/DOI.
