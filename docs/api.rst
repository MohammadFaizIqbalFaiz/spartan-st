API reference
=============

This page summarizes the main public Spartan API. For full parameter explanations,
see the user-guide pages for spatial domain identification, SVG discovery, and
alpha-selection.

Main workflows
--------------

.. py:function:: spartan.tl.spartan_spatial_domains(adata, spatial_coord="grid", spatial_neighs=6, spatial_rings=2, spatial_neighborhood="knn", total_pca_comps=50, pca_comps_extract=30, gene_coord="generic", gene_neighs=15, alpha=0.80, beta1=0.10, beta2=0.40, resolution=1.0, seed=1, key_added="spartan_domains", copy=False)

   Run the complete Spartan spatial domain workflow.

   This function builds the spatial graph, Local Spatial Activation graph,
   gene-expression connectivity graph, aggregated Spartan graph, and then applies
   Leiden clustering. Results are stored in ``adata.obs[key_added]`` and graph
   outputs are stored in ``adata.obsp``.

   See :doc:`spatial_domains` for detailed parameter descriptions.

.. py:function:: spartan.tl.spartan_svg(adata, lsa_graph, layer="log1pX", n_permutations=1000, n_cores=8, use_X_if_missing=True, alpha_svg=0.05, chunk_size=200, seed=1, key_added="spartan_svg", copy=False, dtype=np.float32, prefer_backend="threads", two_stage=True, n_permutations_stage1=100, top_k_refine=3000)

   Run SAQ-based spatially variable gene discovery.

   This function computes Spatial Activation Quotient scores, estimates
   permutation-based p-values, applies Benjamini-Hochberg FDR correction, and
   stores SVG results in ``adata.var``.

   See :doc:`svg` for detailed parameter descriptions.

Preprocessing
-------------

.. py:function:: spartan.tl.pre_process_imaging(adata)

   Preprocess imaging-based spatial transcriptomics datasets such as MERFISH,
   STARmap, osmFISH, BaristaSeq, and Vizgen MERFISH.

   The function stores raw counts in ``adata.layers["raw"]`` and normalized
   log-transformed expression in ``adata.layers["log1pX"]``.

.. py:function:: spartan.tl.pre_process_sequencing(adata, min_counts=500, min_cells=5, jitter=0.4, n_top_genes=3000, seed=1)

   Preprocess sequencing-based spatial transcriptomics datasets such as
   Stereo-seq, Visium, and Visium HD.

   The workflow performs quality control, cell/gene filtering, normalization,
   log transformation, and highly variable gene selection.

Graph construction and clustering
---------------------------------

.. py:function:: spartan.tl.spartan_build_graphs(adata, spatial_coord="grid", spatial_neighs=6, spatial_rings=2, spatial_neighborhood="knn", total_pca_comps=50, pca_comps_extract=30, gene_coord="generic", gene_neighs=15, seed=1, copy=False)

   Build Spartan graph components without running Leiden clustering.

   Stores ``spartan_spatial_graph``, ``spartan_spatial_weights``,
   ``spartan_lsa_graph``, and ``spartan_gene_graph`` in ``adata.obsp``.

.. py:function:: spartan.tl.perform_leiden_clustering(adata, lsa_graph, spatial_graph, gene_graph, alpha, beta1, beta2, resolution, seed, key_added)

   Run Leiden clustering on the aggregated Spartan graph.

   The aggregated graph is stored in ``adata.obsp["spartan_joint_graph"]`` and
   domain labels are stored in ``adata.obs[key_added]``.

Alpha-selection and benchmarking utilities
------------------------------------------

.. py:function:: spartan.tl.initiate_alpha_selection(adata, lower_alpha, upper_alpha, lower_resolution, upper_resolution, step_alpha, step_resolution, lower_nlsas, upper_nlsas, n_jobs, config, seed, use_nLSAS=False)

   Run benchmarking-oriented alpha and resolution operating-regime analysis.

   This workflow assumes ground-truth labels are available for evaluation.
   Ground-truth labels are used for metric calculation, not for Spartan graph
   construction or clustering.

.. py:function:: spartan.tl.consensus_alpha(summary, top_k=50)

   Compute consensus alpha values across multiple samples from alpha-selection
   summaries.

.. py:function:: spartan.tl.compute_nLSAS(lisa_vector, labels)

   Compute the normalized Local Spatial Activation Score for a clustering
   partition.

Plotting API
------------

.. py:function:: spartan.pl.spatial_domains(adata, color="spartan_domains", ax=None, title=None, **kwargs)

   Plot Spartan spatial domains using ``scanpy.pl.spatial``.

.. py:function:: spartan.pl.local_spatial_activation(adata, color="spartan_local_spatial_activation", ax=None, title="Local Spatial Activation", cmap="Reds", **kwargs)

   Plot per-cell or per-spot Local Spatial Activation values.

.. py:function:: spartan.pl.svg_table(adata, n=20, score_key="spartan_saq")

   Return the top ``n`` genes ranked by SAQ score.

.. py:function:: spartan.pl.alpha_selection_summary(summary_df, metric="median_NMI", ax=None)

   Plot alpha-selection summary results.

Internal helper functions
-------------------------

The following functions are exposed for transparency but are mainly intended for
advanced users or developers.

.. py:function:: spartan.tl.build_spatial_graph(adata, crd_type, neighs, rings, neighborhood)

   Build the spatial neighborhood graph and row-normalized spatial weight matrix.

.. py:function:: spartan.tl.perform_pca(adata, comps, comps_extract, seed)

   Perform PCA and return the selected PCA feature matrix.

.. py:function:: spartan.tl.build_lsa_graph(adata, expr_values, spatial_matrix)

   Build the Local Spatial Activation graph.

.. py:function:: spartan.tl.build_gene_expression_connectivity_graph(adata, crd_type, neighs)

   Build the gene-expression connectivity graph in PCA space.

.. py:function:: spartan.tl.to_igraph(adjacency)

   Convert a sparse adjacency matrix into an ``igraph`` graph.

.. py:function:: spartan.tl.normalize(mat)

   Row-normalize a sparse matrix.
