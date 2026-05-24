# Overview

Spartan integrates spatial topology, gene expression connectivity, and Local Spatial Activation (LSA) into an aggregated graph for unsupervised spatial domain detection. It is designed for multiple spatial transcriptomics technologies, including MERFISH, Stereo-seq, and Visium HD.

Spartan constructs three complementary graphs:

- **Spatial graph (`S`)**: physical neighborhood topology.
- **Gene expression connectivity graph (`G`)**: transcriptomic similarity in reduced expression space.
- **Local Spatial Activation graph (`L`)**: neighborhood-conditioned transcriptional activation/deviation structure.

These graphs are combined into an aggregated graph:

```{math}
J = (\alpha-\beta_1)L + (1-\alpha)G + (\alpha-\beta_2)S
```

Leiden clustering is then applied to the aggregated graph to identify spatial domains. Spartan also provides the Spatial Activation Quotient (SAQ) workflow for identifying spatially variable genes using the LSA graph.
