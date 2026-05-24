API reference
=============

This page documents the public Spartan API used by most users. The recommended
user-facing interface is ``spartan.tl`` for analysis functions and ``spartan.pl``
for plotting helpers.

Main workflows
--------------

.. autofunction:: spartan.tl.spartan_spatial_domains

.. autofunction:: spartan.tl.spartan_svg

Preprocessing
-------------

.. autofunction:: spartan.tl.pre_process_imaging

.. autofunction:: spartan.tl.pre_process_sequencing

Graph construction and clustering
---------------------------------

.. autofunction:: spartan.tl.spartan_build_graphs

.. autofunction:: spartan.tl.build_spatial_graph

.. autofunction:: spartan.tl.perform_pca

.. autofunction:: spartan.tl.build_lsa_graph

.. autofunction:: spartan.tl.build_gene_expression_connectivity_graph

.. autofunction:: spartan.tl.perform_leiden_clustering

Alpha-selection and benchmarking utilities
------------------------------------------

.. autofunction:: spartan.tl.initiate_alpha_selection

.. autofunction:: spartan.tl.consensus_alpha

.. autofunction:: spartan.tl.build_graph_cache

.. autofunction:: spartan.tl.mix_weights

.. autofunction:: spartan.tl.compute_nLSAS

.. autofunction:: spartan.tl.prune_results_nLSAS

.. autofunction:: spartan.tl.summarize_all_alphas_nLSAS

Plotting API
------------

.. autofunction:: spartan.pl.spatial_domains

.. autofunction:: spartan.pl.local_spatial_activation

.. autofunction:: spartan.pl.svg_table

.. autofunction:: spartan.pl.alpha_selection_summary

Internal helper functions
-------------------------

These functions are exposed through ``spartan.tl`` for transparency, but most
users do not need to call them directly.

.. autofunction:: spartan.tl.normalize

.. autofunction:: spartan.tl.vec_norm

.. autofunction:: spartan.tl.local_spatial_activation

.. autofunction:: spartan.tl.sparse_to_coord_weight_dict

.. autofunction:: spartan.tl.to_igraph
