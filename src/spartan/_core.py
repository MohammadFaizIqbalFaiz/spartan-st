"""Core SPARTAN methods for spatial domain identification and SVG discovery."""

from __future__ import annotations

import scanpy as sc
from anndata import AnnData
import matplotlib.pyplot as plt
import pandas as pd
import squidpy as sq
import numpy as np
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
from sklearn.metrics import homogeneity_score,completeness_score
import scipy.sparse as sp
from joblib import Parallel, delayed
import igraph as ig
from collections import defaultdict
import leidenalg
from scipy.sparse import coo_matrix
from scipy.sparse import csr_matrix
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
import multiprocessing
from typing import Optional, Tuple, Union


def pre_process_imaging(adata:AnnData):
    """
    Pre-processing function to perform library-size normalization, log-transformation
    on AnnData object, for imaging-based datasets: MERFISH, STARmap, STARmap*, BaristaSeq, OsmFISH, Vizgen MERFISH. 
    The raw count matrix is stored in adata.layers["raw"].
    The final library-size normalized and log-transformed count matrix is stored in adata.layers['log1pX'].

    Returns pre-processed AnnData object.
    """

    adata = adata.copy()
    adata.layers['raw'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['log1pX'] = adata.X.copy()
    adata = adata[adata.obs['ground_truth'].notna()].copy()

    return adata

def pre_process_sequencing(adata:AnnData, min_counts:Optional[int]=500, min_cells:Optional[int]=5, 
                           jitter:Optional[float]=0.4,n_top_genes:Optional[int] = 3000, seed:int=1):
    
    """
    Pre-processing function to compute quality control metrics, and to perform gene and cell filtering,
    library-size normalization, log-transformation on AnnData object, for sequencing-based datasets: 
    Stereo-seq, 10xVisium. 
    The raw count matrix is stored in adata.layers["raw"].
    The final library-size normalized and log-transformed count matrix is stored in adata.layers['log1pX'].

    Parameters:
    
    min_counts: minimum number of counts required for a cell to pass filtering (as in scanpy.pp.filter_cells)
    min_cells: minimum number of cells expressed required for a gene to pass filtering (as in scanpy.pp.filter_genes)
    jitter: jitter value to the stripplot (as in scanpy.pl.violin)
    n_top_genes: Number of highly variable genes to keep, we use the seurat flavor here (as in scanpy.pp.highly_variable_genes)
    seed: The value of seed
    Returns pre-processed AnnData object
    """
    
    adata = adata.copy()
    
    np.random.seed(seed)
    n_genes = adata.shape[1]
    top_n = [n for n in [10, 20, 50, 100] if n < n_genes]

    sc.pp.calculate_qc_metrics(adata, percent_top=top_n, inplace=True)

    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=jitter)

    zero_count_cells = np.where(adata.X.sum(axis=1) == 0)[0]
    print(f"Number of zero-count cells: {len(zero_count_cells)}")

    adata = adata[adata.X.sum(axis=1) > 0].copy()

    # Normalization
    adata.layers['raw'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['log1pX'] = adata.X.copy()
    sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=n_top_genes)

        
    return adata

def normalize(mat): 
   """

   The function performs row normalization of the csr matrix, mat:
   row_means: 1D vector of row means of mat.
   inv_rm: A sprase diagonal matrix where the diagonal elements are the row means.

   Returns a row-normalized sparse matrix mat

   """
   row_means = np.array(mat.sum(axis=1)).flatten()
   row_means[row_means==0]=1
   inv_rm = sp.diags(1.0/row_means)
   mat = inv_rm.dot(mat)

   return mat

def local_spatial_activation(adata:AnnData,mtr):

  """
  The function computes the Local Spatial Activation (LSA) value of each spot/cell:
  It takes the LSA matrix, mtr, and computes the row-sum.

  Parameters:

  adata: an AnnData object.
  mtr: The LSA matrix

  Returns the modified adata with added obs column pertaining to LSA value of each spot/cell

  """

  local_spatial_activation = np.array(mtr.sum(axis=1)).flatten()
  adata.obs['spartan_local_spatial_activation'] = local_spatial_activation

  return 

def vec_norm(v1,v2,flag:str):
    """
    The function computes the L2-norm of the difference of two vectors v1 and v2
    """
    if(v1.shape!=v2.shape):
        raise ValueError("vectors don't align")
    else:
        vt = v1 - v2
        #vt = np.maximum(vt,0.0)
   
    if flag == 'N':
        return np.linalg.norm(vt)
    elif flag == 'S':
        return np.linalg.norm(vt,axis=1)

def build_spatial_graph(adata:AnnData,crd_type:str,neighs:int,rings:int,neighborhood:str):
    
    """
    The functions builds a spatial graph that encodes spatial relationships between two spots/cells. 
    It returns a adjacency matrix S as spat_grh, and a row-normalized spatial weights matrix W^r as spatial_W

    Input: 
    
    adata: AnnData object.
    crd_type: Type of coordinate system. "grid" for sequencing-based datasets such as 10x Visium/Stereo-seq/Visium HD.
              "generic" for imaging-based datasets such as MERFISH.
    neighs: The number of neighboring spots/cells inside a neighborhood.
    rings: The number of rings of neighbors for grid data.
    neighborhood: Type of spatial graph. "knn" for KNN graph. "delaunay" for Delaunay graph.

    Output:

    adata: Updated AnnData object.
    spat_grh: The adjacency matrix of the spatial graph.
    spatial_W: Row-normalized spatial weights matrix.

    """
    
    #Creates spatial neighborhood graphs based on spatial distances (geographical locations)
    if neighborhood == "knn":
        sq.gr.spatial_neighbors(adata,coord_type=crd_type,n_neighs=neighs,n_rings=rings,key_added = 'spatial_graph')
    else:    
        sq.gr.spatial_neighbors(adata, delaunay=True, coord_type=crd_type, key_added = 'spatial_graph')
        

    spat_grh = adata.obsp['spatial_graph_connectivities']

    #row-standarised spatial weight matrix
    spatial_W = adata.obsp['spatial_graph_distances']#Gets spatial distance matrix
    epsilon = 1e-3
    spatial_W.data = 1.0/(spatial_W.data + epsilon)
    spatial_W = normalize(spatial_W)

    adata.obsp['spartan_spatial_distances'] = spatial_W #Assigns row-standarized spatial weight matrix

    return adata,spat_grh,spatial_W

def perform_pca(adata:AnnData,comps:int,comps_extract:int,seed:int):

    """
    The function performs PCA on the filtred-expression matrix X and transforms it into feature matrix Y.

    Input:

    adata: AnnData object.
    comps: Number of principal components to compute.
    comps_extract: Number of principal components to keep.
    seed: Seed value to ensure reproducibilty
    
    Output:

    adata: Updated AnnData object
    PCA: A feature matrix where the number of columns is equal to the comps_extract.

    """
    
    #Perform PCA

    sc.pp.scale(adata)#Centering of the data
    sc.pp.pca(adata, n_comps=comps, svd_solver='arpack',random_state = seed) #Performs PCA

    PCA = adata.obsm['X_pca'][:, :comps_extract]#Extracts PCA values

    return adata,PCA

def build_lsa_graph(adata:AnnData,expr_values,spatial_matrix):

   """
   The sparse edge-list verison of the LSA graph formulation that builds the LSA graph and constructs the weighted adjacency matrix L. 
   Each entry l_ij is the LSA weight of the spatial edge i -> j. 

   Input:

   adata: AnnData object.
   expr_values: The transformed featured matrix where rows denote spots and columns
                denote the number of PCs/features.
   spatial_matrix: Row-normalized spatial weights matrix.

   Output:

   G: The weighted adjacency matrix of the LSA graph. 

   """

   
   G = spatial_matrix.tocoo()#Converts it into COO format

   local_means = np.zeros_like(expr_values)

   local_vars = np.zeros(adata.n_obs)

   for i in range(adata.n_obs):
        neighbors = G.row == i
        indices_i = G.col[neighbors]
        all_indices = np.append(indices_i,i)
        if np.any(neighbors):
          local_means[i] = np.mean(expr_values[all_indices],axis=0)
          vc = np.var(expr_values[all_indices],axis=0)
          local_vars[i]= np.linalg.norm(vc)
    
    
    
   row = G.row #The array of source nodes
   col = G.col #The array of target nodes
   data = G.data #The array of spatial weights associated with edges (source -> target)

   xi = expr_values[row] #The matrix where each row is the feature vector of the source node
   xj= expr_values[col]  #The matrix where each row is the feature vector of the target node

   local_mean_i = local_means[row] #The matrix where each row is the local neighborhood mean vector of the source node
   local_var_i = local_vars[row] #The matrix where each row is the local neighborhood variance vector of the source node

   diffI = vec_norm(xi,local_mean_i,flag='S') #The diffI matrix where each row is the deviation vector of the source node 
                                              #from its local neighborhood
   diffI = np.nan_to_num(diffI, nan=0.0) #Converts NaN values to 0.0
   diffJ = vec_norm(xj,local_mean_i,flag='S') #The diffJ matrix where each row is the deviation vector of the target node 
                                              #from the local neighborhood of the source node
   diffJ = np.nan_to_num(diffJ, nan=0.0) #Converts NaN values to 0.0

   attribute_deviation = diffI * diffJ #The matrix where each entry [i,j] is the unormalized joint-deviation weight of the spatial edge i->j
   attribute_deviation = attribute_deviation / (local_var_i + np.finfo(float).eps) #The matrix where each entry [i,j] is the normalized joint-deviation weight 
                                                                                   #of the spatial edge i->j
   adjusted_weight = data * attribute_deviation  #The matrix where each entry [i,j] is the normalized LSA weight of the spatial edge i->j
   G = coo_matrix((adjusted_weight, (row, col)), shape=G.shape) #The LSA matrix in coordinate form.

   G = G.tocsr() #Converts the matrix back to sparse format

   return G

def build_gene_expression_connectivity_graph(adata:AnnData,crd_type:str,neighs:int):

    """
    The function builds a gene expression connectivity graph.

    Input:

    adata: AnnData object.
    crd_type: The type of coordinate system (grid or generic).
    neighs: The number of neighbors to consider to build the KNN graph.

    Output:

    gene_graph: The adjacency matrix of the gene expression connectivity graph.

    """
    
    #Create neighborhood graphs based on gene expression similarity in PCA coordinate space
    sq.gr.spatial_neighbors(adata, coord_type=crd_type, n_neighs=neighs, key_added = 'gene_graph')
    gene_graph = adata.obsp['gene_graph_connectivities']

    return gene_graph

def spartan_build_graphs(adata:AnnData,spatial_coord:str="grid",spatial_neighs:Optional[int]=6,spatial_rings:Optional[int]=2,
                            spatial_neighborhood:str="knn",
                            total_pca_comps:int=50,pca_comps_extract:int=30,gene_coord:str="generic",gene_neighs:int=15,seed:int=1,
                            copy: bool=False):

    """
    The function is a wrapper that calls the functions to build the spatial graph, Local Spatial Activation graph,
    gene expression connectivity graph.

    Input:

    adata: AnnData object.
    spatial_coord: The type of coordinate system for the spatial graph (default: grid).
    spatial_neighs: The number of neighbors to consider to build the spatial graph (default: 6).
    spatial_rings: The number of rings to consider to build the spatial graph(only required if spatial_coord = "grid", default:2).
    spatial_neighborhood: The type of neighborhood to consider to build the spatial graph (KNN or Delaunay, default:"knn").
    total_pca_comps: The total number of PCs to compute for PCA (default:50).
    pca_comps_extract: The number of PCs to keep (default:30).
    gene_coord: The type of coordinate system for the gene expression connectivity graph (default: generic).
    gene_neighs: The number of neighbors to consider to build the gene expression connectivity graph (default: 15).
    seed: The seed value to ensure reproducibilty (default: 1)
    copy: Whether to create a copy of the AnnData object. (Default: False)

    Output:

    adata: If copy = True. Returns the updated AnnData object:
            adata.obsp['spartan_spatial_graph'], The obsp column that keeps the adjacency matrix of the spatial graph.
            adata.obsp['spartan_spatial_weights'], The obsp column that keeps the row-normalized spatial weights matrix.
            adata.obsp['spartan_lsa_graph'], The obsp column that keeps the weighted adjacency matrix of the LSA graph.
            adata.obsp['spartan_gene_graph'], The obsp column that keeps the adjacency matrix of the gene expression connectivity graph.

            If copy = False, the all the updates are done on the original AnnData object.

    """
           

    if copy:
        adata=adata.copy()

    #Builds the spatial graph    
    
    adata,spatial_graph,spatial_W = build_spatial_graph(adata,crd_type=spatial_coord,neighs=spatial_neighs,rings=spatial_rings,
                                                       neighborhood=spatial_neighborhood)
    adata.obsp['spartan_spatial_graph'] = spatial_graph
    adata.obsp['spartan_spatial_weights'] = spatial_W

    #Performs PCA
    
    adata,PCA = perform_pca(adata,comps=total_pca_comps,comps_extract=pca_comps_extract,seed=seed)

    expression_values = PCA.toarray() if hasattr(PCA, 'toarray') else PCA#Creates a dense array

    #Builds the LSA graph

    lsa_graph = build_lsa_graph(adata,expr_values=expression_values,spatial_matrix=spatial_W)

    adata.obsp['spartan_lsa_graph'] = lsa_graph
    
    #construct gene expression connectivity graph
    spac = adata.obsm['spatial']
    adata.obsm['spatial'] = PCA #Uses PCA coordinates to create gene expression similarity based neighborhood graphs

    try:
      gene_graph = build_gene_expression_connectivity_graph(adata,crd_type=gene_coord,neighs=gene_neighs)
      adata.obsp['spartan_gene_graph']=gene_graph
    finally:
      adata.obsm["spatial"]=spac
        
    if copy:
        return adata
    return None

#Convert sparse matrix to coo format

def sparse_to_coord_weight_dict(mat):

    """
    The function transforms the matrix in CSR format to COO format and creates a dictionary 
    where each element is of the form {(source,target):weight}. 

    Input:
    mat: CSR matrix

    Output:
    d: Returns the matrix mat in dictionary format.

    """
     

    coo = mat.tocoo()

    d = {}

    for i,j,w in zip(coo.row,coo.col,coo.data):
        d[(int(i),int(j))] = float(w)

    return d

def build_graph_cache(lsa_graph=None,gene_graph=None,spatial_graph=None):

    """
    The function creates a template of the aggregated graph that aggregates spatial graph, LSA graph, and 
    the gene expression connectivity graph. Note on asymmetric Spartan graphs:
    The aggregated graph J may be asymmetric because the LSA graph is neighborhood-conditioned, 
    so J_ij and J_ji can encode different local activation relationships. For spatial domain clustering, 
    Spartan intentionally converts J to an undirected igraph without simplifying reciprocal edges. 
    Thus reciprocal entries are retained as parallel edges, allowing Leiden to use their 
    combined bidirectional affinity: J_eff(i,j) = J_ij + J_ji.
    Do not symmetrize or simplify this graph unless intentionally changing the validated Spartan clustering behavior.


    Input:

    lsa_graph: The weighted adjacency matrix of the LSA graph in CSR format.
    gene_graph: The adjacency matrix of the gene expression connectivity graph in CSR format.
    spatial_graph: The adjacency matrix of the spatial graph in CSR format.

    Output:

    Returns a dictionary of object:
    g: The template of the aggregated graph (igraph object) where vertices denote spots/cells.
       The edges are the set of union of edges belonging to spatial graph, LSA graph, and gene expression connectivity graph.
    coords: The set of union of edges belonging to spatial graph, LSA graph, and gene expression connectivity graph.
    w_lsa: The array of edge weights associated with the LSA graph.
    w_gene: The array of edge weights associated with the gene expression connectivity graph.
    w_lsa: The array of edge weights associated with the spatial graph.

    """

    layer_dicts={}

    if lsa_graph is not None:
        layer_dicts["lsa"] = sparse_to_coord_weight_dict(lsa_graph)
        n_nodes = lsa_graph.shape[0]

    if gene_graph is not None:
        layer_dicts["gene"] = sparse_to_coord_weight_dict(gene_graph)
        n_nodes = gene_graph.shape[0]

    if spatial_graph is not None:
        layer_dicts["spatial"] = sparse_to_coord_weight_dict(spatial_graph)
        n_nodes = spatial_graph.shape[0]


    all_coords = sorted(set().union(*[d.keys() for d in layer_dicts.values()]))

    sources = np.array([i for i,j in all_coords], dtype=np.int32)
    targets = np.array([j for i,j in all_coords], dtype=np.int32)

    g = ig.Graph(directed=False)
    g.add_vertices(n_nodes)
    g.add_edges(list(zip(sources,targets)))

    w_lsa = np.array([layer_dicts.get("lsa",{}).get(coord,0.0) for coord in all_coords], dtype=np.float64)
    w_gene = np.array([layer_dicts.get("gene",{}).get(coord,0.0) for coord in all_coords], dtype=np.float64)
    w_spatial = np.array([layer_dicts.get("spatial",{}).get(coord,0.0) for coord in all_coords], dtype=np.float64)

    return {

      "graph": g,
      "coords": all_coords,
      "w_lsa": w_lsa,
      "w_gene": w_gene,
      "w_spatial": w_spatial
    }

def mix_weights(alpha, config, w_lsa=None, w_gene=None, w_spatial=None):

    """
    The function builds one of the three different configuations of the aggregated graph based on the value of config.
    The "lsg" denotes the standard configuration of Spartan that aggregates spatial graph, LSA graph, 
    and gene expression connectivity graph. The "sg" denotes only spatial and gene expression connectivity. 
    The "lg" denotes only LSA and gene expression connectivity. 
    Beta1 and Beta2 are set to default here 0.26 and 0.24 for alpha selection in benchmarking study.

    Input:

    alpha: The primary multiplex weighting parameter. 
    config: The aggregate graph configuration.
    w_lsa: The array of edge weights associated with the LSA graph.
    w_gene: The array of edge weights associated with the gene expression connectivity graph.
    w_spatial: The array of edge weights associated with the spatial graph.

    Output:

    Returns the array of edge weights associated with the aggregated graph.

    """
    
    if config == "lsg":
        return (alpha - 0.26) * w_lsa + (1 - alpha) * w_gene + (alpha - 0.24) * w_spatial
    elif config == "lg":
        return alpha * w_lsa + (1 - alpha) * w_gene
    elif config == "sg":
        return alpha * w_spatial + (1 - alpha) * w_gene
    else:
        raise ValueError("Unknown config")        

def compute_nLSAS(lisa_vector,labels):

    """
    The function computes the normalized Local Spatial Activation Score (nLSAS) of a clustering partition.

    Input:

    lisa_vector: The vector of LSA values of each spot/cell corresponding to the LSA graph.
    labels: The clustering labels of each spot/cell.

    Output:

    Returns the average normalized activation across all clusters of the given partition.

    """

    labels = np.array(labels)
    unique_clusters = np.unique(labels)
    lisa_scores = {}
    for k in  unique_clusters:
        mask = [labels==k]
        mask = np.array(mask,ndmin=1)
        len_cluster = len(mask==1)
        lisa_vec = np.dot(lisa_vector,mask.T) 
        lisa_scores[k] = np.mean(lisa_vec)

    max_score = max(lisa_scores.values())
    min_score = min(lisa_scores.values())

  
    normalized_lisa = {k: ((v - min_score) / (max_score - min_score + 1e-10))   for k, v in lisa_scores.items()}
    
    return pd.Series(normalized_lisa, name="normalized_lisa_score").sort_index() 


def process_alpha_light(
    alpha,
    graph_template,
    resolutions,
    config,
    true_labels,
    seed,
    w_lsa=None,
    w_gene=None,
    w_spatial=None,
    nLSAS_values=None,
    use_nLSAS=False,
):

    """
    For a given alpha, the function computes the evaluation metrics NMI, HOM, HOM, COM, nLSAS of 
    partitions (clustering assignments from the Leiden clustering) for every resolution across the resolution grid.

    Input:

    alpha: The primary multiplex weighting parameter.
    graph_template: The igraph object of the aggregated graph.
    resolutions: The resolution grid.
    config: The configuration of the aggregated graph (lsg/sg/lg)
    true_labels: ground truth labels.
    seed: the seed value.
    w_lsa: the array of edge weights associated with the LSA graph.
    w_gene: the array of edge weights associated with the gene expression connectivity graph.
    w_spatial: the array of edge weights associated with the spatial graph.
    nLSAS_values: the vector of LSA values of each spot/cell (intialization).
    use_nLSAS: Whether to compute the nLSAS of a partition.

    Output:

    Returns a list, alpha_results which contains for every resolution{
    alpha: the alpha parameter.
    resolution: the value of the resolution.
    ClusterN: The size of the partition.
    NMI: Normalized Mutual Information (NMI) score of the partition.
    HOM: Homogeneity (HOM) score of the partition.
    COM: Completeness (COM) score of the partition.
    nLSAS: normalized Local Spatial Activation Score (nLSAS) of the partition.

    """
    
    g = graph_template.copy()

    weights = mix_weights(
        alpha=alpha,
        config=config,
        w_lsa=w_lsa,
        w_gene=w_gene,
        w_spatial=w_spatial,
    )
    g.es["weight"] = weights

    alpha_results = []

    for res in resolutions:
        partition = leidenalg.find_partition(
            g,
            leidenalg.RBConfigurationVertexPartition,
            weights=g.es["weight"],
            n_iterations=-1,
            resolution_parameter=float(res),
            seed=seed,
        )

        prt = np.asarray(partition.membership)
        labels = np.unique(prt)

        row = {
            "alpha": round(float(alpha), 3),
            "resolution": round(float(res), 3),
            "ClusterN": len(labels),
            "NMI": normalized_mutual_info_score(true_labels, prt),
            "HOM": homogeneity_score(true_labels, prt),
            "COM": completeness_score(true_labels, prt),
        }

        if use_nLSAS and nLSAS_values is not None:
            nLSAS = compute_nLSAS(nLSAS_values, prt)
            row["nLSAS"] = nLSAS.mean()
        alpha_results.append(row)

    return alpha_results

def run_grid_search_parallel_alpha_light(
    alpha_grid,
    resolution_grid,
    config,
    cache,
    true_labels,
    seed=0,
    nLSAS_values=None,
    use_nLSAS=False,
    n_jobs=4,
):

    """

    The function is a wrapper for the function process_alpha_light(). It intializes the grids for alpha and resolution parameters
    and executes the function process_alpha_light() in parallel across multiple threads for each alpha and the resolution grid.

    Input:

    alpha_grid: The grid for the alpha parameter.
    resolution_grid: The grid for the resolution parameter.
    config: The configuration of the aggregated graph.
    cache: The dictionary that contains the template for the aggregated, spatial, LSA, and gene expression connectivity graphs.
    true_labels: ground truth labels.
    seed: the seed value.
    nLSAS_values: the vector of LSA values of each spot/cell (intialization).
    use_nLSAS: Whether to compute the nLSAS of a partition.
    n_jobs: Number of threads to use.

    Output:

    Returns a dataframe that conatins the results of the function process_alpha_light(), for every alpha across the alpha grid.

    """
    alpha_grid = np.round(np.asarray(alpha_grid, dtype=float), 3)
    resolution_grid = np.round(np.asarray(resolution_grid, dtype=float), 3)

    parallel_out = Parallel(n_jobs=n_jobs)(
        delayed(process_alpha_light)(
            alpha=alpha,
            graph_template=cache["graph"],
            resolutions=resolution_grid,
            config=config,
            true_labels=true_labels,
            seed=seed,
            w_lsa=cache["w_lsa"],
            w_gene=cache["w_gene"],
            w_spatial=cache["w_spatial"],
            nLSAS_values=nLSAS_values,
            use_nLSAS=use_nLSAS,
        )
        for alpha in alpha_grid
    )

    results = [row for alpha_rows in parallel_out for row in alpha_rows]
    df_results = pd.DataFrame(results)

    return df_results

def prune_results_nLSAS(
    adata,
    df,
    alpha,
    lower_nLSAS_percentile=40,
    upper_nLSAS_percentile=70,
):
    """
    For one alpha:
    1. keep only rows with ClusterN == true domain count
    2. apply nLSAS percentile filtering
    3. return the pruned dataframe
    """
    num_domains = len(np.unique(adata.obs["ground_truth"]))

    # keep rows for this alpha
    df_alpha = df[np.isclose(df["alpha"], alpha)].copy()

    if df_alpha.empty:
        return pd.DataFrame()

    # keep only rows with correct cluster count
    df_alpha = df_alpha[df_alpha["ClusterN"] == num_domains].copy()

    if df_alpha.empty:
        return pd.DataFrame()

    # drop rows where nLSAS is missing
    df_alpha = df_alpha.dropna(subset=["nLSAS"]).copy()

    if df_alpha.empty:
        return pd.DataFrame()

    # compute percentile thresholds within this alpha
    lower_thr = np.percentile(df_alpha["nLSAS"], lower_nLSAS_percentile)
    upper_thr = np.percentile(df_alpha["nLSAS"], upper_nLSAS_percentile)

    # retain only rows in the desired nLSAS band
    df_alpha = df_alpha[
        (df_alpha["nLSAS"] >= lower_thr) &
        (df_alpha["nLSAS"] <= upper_thr)
    ].copy()

    return df_alpha

def summarize_all_alphas_nLSAS(
    adata,
    df_results,
    lower_nLSAS_percentile=40,
    upper_nLSAS_percentile=70,
):
    """
    For each alpha:
    - prune by ClusterN and nLSAS band
    - compute median NMI/HOM/COM/nLSAS
    - record best config by NMI
    """
    results_summary = []

    for alpha in sorted(df_results["alpha"].dropna().unique()):
        df_pruned = prune_results_nLSAS(
            adata=adata,
            df=df_results,
            alpha=alpha,
            lower_nLSAS_percentile=lower_nLSAS_percentile,
            upper_nLSAS_percentile=upper_nLSAS_percentile,
        )

        if df_pruned.empty:
            continue

        # median metrics after pruning
        median_stats = df_pruned[["NMI", "HOM", "COM", "nLSAS"]].median(axis=0)

        # best config within this alpha, using NMI
        best_config = df_pruned.loc[df_pruned["NMI"].idxmax()]

        results_summary.append({
            "alpha": alpha,
            "median_NMI": median_stats["NMI"],
            "median_HOM": median_stats["HOM"],
            "median_COM": median_stats["COM"],
            "median_nLSAS": median_stats["nLSAS"],
            "best_NMI": best_config["NMI"],
            "best_HOM": best_config["HOM"],
            "best_COM": best_config["COM"],
            "best_nLSAS": best_config["nLSAS"],
            "best_resolution": best_config["resolution"],
            "best_config_index": best_config.name,
            "n_retained": len(df_pruned),
        })

    summary_df = pd.DataFrame(results_summary)

    if not summary_df.empty:
        summary_df = summary_df.sort_values(
            by="median_NMI",
            ascending=False
        ).reset_index(drop=True)

    return summary_df

def initiate_alpha_selection(adata:AnnData,lower_alpha:float,upper_alpha:float,
                            lower_resolution:float,upper_resolution:float,step_alpha:float,step_resolution:float,lower_nlsas:int,upper_nlsas:int,
                            n_jobs:int,config:str,seed:int,use_nLSAS:bool = False):

    """
    The function intiates the alpha selection strategy. alpha_grid = [lower_alpha,upper_alpha].
    resolution_grid = [lower_resolution,upper_resolution]. 
    The function assumes the ground truth exists as it is the case in benchmarking study.

    Input:

    adata: AnnData object.
    lower_alpha: The lower bound of the alpha parameter.
    upper_alpha: The upper bound of the alpha parameter.
    lower_resolution: The lower bound of the resolution parameter.
    upper_resolution: The upper bound of the resolution parameter.
    step_alpha: The step size for the alpha grid. 
    step_resolution: The step size for the resolution grid.
    lower_nlsas: The lower bound for the nLSAS metric.
    upper_nlsas: The upper bound for the nLSAS metric.
    n_jobs: The number of threads to consider for parallelization.
    config: The configuration of the aggegated graph (lsg = LSA+Spatial+Gene, lg = LSA+Gene, sg = Spatial+Gene)
    seed: The value of the seed.
    use_nLSAS: Whether to compute the nLSAS metric (True for lsg and lg configurations).

    Output:

    Returns the dataframe that contains the summary of results for the given alpha grid and resolution grid.

    """
    
        
    spatial_graph = adata.obsp["spartan_spatial_graph"]
    lsa_graph = adata.obsp["spartan_lsa_graph"]
    gene_graph = adata.obsp["spartan_gene_graph"]
    
    true_labels = np.asarray(adata.obs["ground_truth"])
    local_spatial_activation(adata,lsa_graph)
    nLSAS_values = np.asarray(adata.obs["spartan_local_spatial_activation"])

    alpha_grid = np.arange(lower_alpha, upper_alpha, step_alpha)
    resolution_grid = np.arange(lower_resolution, upper_resolution, step_resolution)

    cache = build_graph_cache(
            lsa_graph=lsa_graph,
            gene_graph=gene_graph,
            spatial_graph=spatial_graph,
            )

    df_results = run_grid_search_parallel_alpha_light(
                 alpha_grid=alpha_grid,
                 resolution_grid=resolution_grid,
                 config=config,
                 cache=cache,
                 true_labels=true_labels,
                 seed=seed,
                 nLSAS_values=nLSAS_values,
                 use_nLSAS=use_nLSAS,
                 n_jobs=n_jobs,
                 )

    #print(df_results.head())
    
    summary_df = summarize_all_alphas_nLSAS(
    adata=adata,
    df_results=df_results,
    lower_nLSAS_percentile=lower_nlsas,
    upper_nLSAS_percentile=upper_nlsas,
    )

    #print(summary_df.head())

    return summary_df

    

def consensus_alpha(summary: list[pd.DataFrame], top_k: int = 50):

    """
    The function computes the top k alphas that perform well in terms of median NMI across 
    all the samples of a given spatially resolved trancriptomics (SRT).

    Input:

    summary: The list of datframes that conatins the results for all samples (for a given alpha grid and a resolution grid)
    top-k: The number of best performing alphas to consider (default:50)

    Output:

    Returns the dataframe that contains the results of the top k alphas (sorted by the mean of the median NMI) across all samples.

    """

    if len(summary) == 0:
        return pd.DataFrame()

    top_alpha_sets = []

    for df in summary:

        top_alphas = df.sort_values(by='median_NMI', ascending=False).head(top_k)["alpha"]
        top_alpha_sets.append(set(top_alphas))

    common_alphas = set.intersection(*top_alpha_sets)

    if len(common_alphas) > 0:
        selected_alphas = common_alphas
    else:
        selected_alphas = set.union(*top_alpha_sets)

    renamed_dfs = []

    for i,df in enumerate(summary,start=1):

        sub = df[df["alpha"].isin(selected_alphas)][["alpha","median_NMI"]].copy()
        sub = sub.rename(columns = {"median_NMI": f"median_NMI_s{i}"})
        renamed_dfs.append(sub)


    combined = renamed_dfs[0]
    for df in renamed_dfs[1:]:

         combined = combined.merge(df, on="alpha", how="outer")

    nmi_cols = [col for col in combined.columns if col.startswith("median_NMI_s")]
    combined["mean_median_NMI"] = combined[nmi_cols].mean(axis=1)

    combined = combined.sort_values(by="mean_median_NMI", ascending=False).reset_index(drop=True)

    return combined

def to_igraph(adjacency):

    """
    The function creates an igraph object of the aggregated graph. Note on asymmetric Spartan graphs:
    The aggregated graph J may be asymmetric because the LSA graph is neighborhood-conditioned, 
    so J_ij and J_ji can encode different local activation relationships. For spatial domain clustering, 
    Spartan intentionally converts J to an undirected igraph without simplifying reciprocal edges. 
    Thus reciprocal entries are retained as parallel edges, allowing Leiden to use their 
    combined bidirectional affinity: J_eff(i,j) = J_ij + J_ji.
    Do not symmetrize or simplify this graph unless intentionally changing the validated Spartan clustering behavior.

    Input:

    adjacency: The weighted adjacency matrix of the graph in sparse format.

    Output:

    Returns the igraph object g.
    
    """
    
    sources, targets = adjacency.nonzero()
    weights = adjacency[sources, targets].A1  # Convert sparse matrix to flat array
    g = ig.Graph(directed=False)
    g.add_vertices(adjacency.shape[0])
    g.add_edges(zip(sources, targets))
    g.es['weight'] = weights
    return g

def perform_leiden_clustering(adata:AnnData,lsa_graph:csr_matrix,spatial_graph:csr_matrix,gene_graph:csr_matrix,
                              alpha:float,beta1:float,beta2:float,resolution:float,seed:int,key_added:str):

    """
    The function executes leiden clustering on the the aggregated graph (joint_graph). 
    It builds an aggregated graph by adding Spatial (spatial_graph), LSA (lsa_graph), and 
    gene expression connectivity (gene_graph) graphs.

    Input:

    adata: AnnData object.
    lsa_graph: The weighted adjacency matrix of the LSA graph.
    spatial_graph: The adjacency matrix of the spatial graph.
    gene_graph: The adjacency matrix of the gene expression connectivity graph.
    alpha: The primary multiplex weighting parameter.
    beta1: The parameter that determines the contribution (alpha-beta1) of the LSA graph.
    beta2: The parameter that determines the contribution (alpha-beta2) of the spatial graph.
    seed: The value of the seed.
    key_added: The key added to the newly created obs columns of the AnnData objects.

    Output:

    Returns the updated AnnData object with newly created obs column adata.obs['spartan_leiden'],
    and adata.obs['key_added'] that contains the partition of the aggregated graph as category type.

    """

    joint_graph = (alpha-beta1) * lsa_graph + (1 - alpha) * gene_graph + (alpha-beta2) * spatial_graph

    adata.obsp['spartan_joint_graph'] = joint_graph
    
    j_graph = to_igraph(joint_graph)

    partition = leidenalg.find_partition(
        j_graph,
        leidenalg.RBConfigurationVertexPartition,
        weights=j_graph.es['weight'],n_iterations=-1,seed=seed,
        resolution_parameter= resolution)

    adata.obs['spartan_leiden'] = [str(c) for c in partition.membership]
    adata.obs[key_added] =  adata.obs['spartan_leiden'].astype('category')

    return adata

def spartan_spatial_domains(
    adata: AnnData,
    spatial_coord: str = "grid",
    spatial_neighs: Optional[int] = 6,
    spatial_rings: Optional[int] = 2,
    spatial_neighborhood: str = "knn",
    total_pca_comps: int = 50,
    pca_comps_extract: int = 30,
    gene_coord: str = "generic",
    gene_neighs: int = 15,
    alpha: float = 0.80,
    beta1: float = 0.1,
    beta2: float = 0.4,
    resolution: float = 1.0,
    seed: int = 1,
    key_added: str = "spartan_domains",
    copy: bool = False,
):

    """
    The function is the wrapper that calls the functions spartan_build_graphs() and perform_leiden_clustering().
    It is the primary function to be used if spatial domains are to be identified using user-defined parameters.
    The beta1 + beta2 must be 0.5.

    Input:

    adata: AnnData object.
    spatial_coord: The type of coordinate system for the spatial graph (default: grid).
    spatial_neighs: The number of neighbors to consider to build the spatial graph (default: 6).
    spatial_rings: The number of rings to consider to build the spatial graph(only required if spatial_coord = "grid", default:2).
    spatial_neighborhood: The type of neighborhood to consider to build the spatial graph (KNN or Delaunay, default:"knn").
    total_pca_comps: The total number of PCs to compute for PCA (default:50).
    pca_comps_extract: The number of PCs to keep (default:30).
    gene_coord: The type of coordinate system for the gene expression connectivity graph (default: generic).
    gene_neighs: The number of neighbors to consider to build the gene expression connectivity graph (default: 15).
    alpha: The primary multiplex weighting parameter (alpha belongs to (0,1)).
    beta1: The parameter that determines the contribution (alpha-beta1) of the LSA graph (beta1 <= alpha).
    beta2: The parameter that determines the contribution (alpha-beta2) of the spatial graph (beta2 <= alpha).
    resolution: The resolution parameter for the leiden clustering.
    seed: The value of the seed.
    key_added: The key added to the newly created obs columns of the AnnData objects.
    copy: Whether to create a copy of the AnnData object. (Default: False)

    Output:

    adata: If copy = True. Returns the updated AnnData object:
            adata.obsp['spartan_spatial_graph'], The obsp column that keeps the adjacency matrix of the spatial graph.
            adata.obsp['spartan_spatial_weights'], The obsp column that keeps the row-normalized spatial weights matrix.
            adata.obsp['spartan_lsa_graph'], The obsp column that keeps the weighted adjacency matrix of the LSA graph.
            adata.obsp['spartan_gene_graph'], The obsp column that keeps the adjacency matrix of the gene expression connectivity graph.
            adata.obs['key_added']: The Spartan's spatial domains.

            If copy = False, all updates are done on the original AnnData object.

    """
    
    if copy:
        adata = adata.copy()

    # Parameter checks
    if not (0 < alpha < 1):
        raise ValueError("alpha must be between 0 and 1.")

    if beta1 < 0 or beta2 < 0:
        raise ValueError("beta1 and beta2 must be non-negative.")

    if beta1 > alpha or beta2 > alpha:
        raise ValueError("beta1 and beta2 must be less than or equal to alpha.")

    if not np.isclose(beta1 + beta2, 0.5):
        raise ValueError("beta1 and beta2 must sum to 0.5.")

    if resolution <= 0:
        raise ValueError("resolution must be positive.")

    if pca_comps_extract > total_pca_comps:
        raise ValueError("pca_comps_extract cannot exceed total_pca_comps.")

    if spatial_neighborhood not in {"knn", "delaunay"}:
        raise ValueError("spatial_neighborhood must be either 'knn' or 'delaunay'.")

    if spatial_coord == "grid" and spatial_rings is None:
        raise ValueError("spatial_rings is required when spatial_coord='grid'.")

    if spatial_neighs is not None and spatial_neighs <= 0:
       raise ValueError("spatial_neighs must be positive.")

    if gene_neighs <= 0:
        raise ValueError("gene_neighs must be positive.")

    # Build graphs
    spartan_build_graphs(
        adata,
        spatial_coord=spatial_coord,
        spatial_neighs=spatial_neighs,
        spatial_rings=spatial_rings,
        spatial_neighborhood=spatial_neighborhood,
        total_pca_comps=total_pca_comps,
        pca_comps_extract=pca_comps_extract,
        gene_coord=gene_coord,
        gene_neighs=gene_neighs,
        seed=seed,
        copy=False,
    )

    #Store the graphs
    lsa_graph = adata.obsp['spartan_lsa_graph']
    spatial_graph = adata.obsp['spartan_spatial_graph']
    gene_graph = adata.obsp['spartan_gene_graph']
    
    # Perform clustering
    perform_leiden_clustering(
        adata,
        lsa_graph=lsa_graph,
        spatial_graph= spatial_graph,
        gene_graph= gene_graph,
        alpha=alpha,
        beta1=beta1,
        beta2=beta2,
        resolution=resolution,
        seed=seed,
        key_added=key_added,
    )

    if copy:
        return adata

    return None

def spartan_svg(
    adata: AnnData,
    lsa_graph: csr_matrix,
    layer: str | None = "log1pX",
    n_permutations: int = 1000,
    n_cores: int = 8,
    use_X_if_missing: bool = True,
    alpha_svg: float = 0.05,
    chunk_size: int = 200,
    seed: int = 1,
    key_added: str = "spartan_svg",
    copy: bool = False,
    dtype=np.float32,
    prefer_backend: str = "threads",   # "threads" or "processes"

    # ----------------------------
    # Two-stage refinement controls
    # ----------------------------
    two_stage: bool = True,
    n_permutations_stage1: int = 100,  # small run for all genes
    top_k_refine: int = 3000,          # refine only these genes
):
    """
    Spartan SAQ/SVG scoring with permutation-based normal approximation.

    Observed score per gene j:
      SAQ_j = (x'_j^T L x'_j) / ||x'_j||^2
    where x'_j is mean-centered expression across spots for gene j.

    P-values:
      We estimate null mean/variance of SAQ_j by permutations (row shuffles),
      then compute z = (obs - mean_null) / std_null and one-sided p = sf(z).

    Two-stage mode (recommended if you want top_k_refine speedups):
      Stage 1: n_permutations_stage1 for all genes -> p-values for all genes
      Stage 2: additional permutations only for top_k_refine genes -> refined p-values
      Final BH/FDR uses final p-values for all genes.
    """
    if copy:
        adata = adata.copy()

    # --------------------------
    # 0) Load X
    # --------------------------
    if layer is not None and layer in adata.layers:
        X = adata.layers[layer]
    elif use_X_if_missing:
        X = adata.X
    else:
        raise KeyError(
            f"Layer {layer!r} not found in adata.layers and use_X_if_missing=False."
        )

    # --------------------------
    # 1) Ensure sparse formats
    # --------------------------
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    else:
        X = X.tocsr()

    if not sp.issparse(lsa_graph):
        lsa_graph = sp.csr_matrix(lsa_graph)
    else:
        lsa_graph = lsa_graph.tocsr()

    N, G = X.shape

    if n_permutations < 1:
        raise ValueError("n_permutations must be >= 1")

    if n_cores < 1:
        n_cores = 1

    # Two-stage sanity
    if two_stage:
        if n_permutations_stage1 < 1:
            raise ValueError("n_permutations_stage1 must be >= 1")
        if n_permutations_stage1 >= n_permutations:
            # no stage2 to do
            two_stage = False
        top_k_refine = int(np.clip(top_k_refine, 0, G))

    # --------------------------
    # 2) Output arrays
    # --------------------------
    all_obs = np.zeros(G, dtype=np.float64)
    all_pvals = np.ones(G, dtype=np.float64)
    stage1_pvals = None
    if two_stage:
        stage1_pvals = np.ones(G, dtype=np.float64)

    # --------------------------
    # 3) Precompute gene means (sparse-safe)
    # --------------------------
    X_mean = np.asarray(X.mean(axis=0)).ravel().astype(np.float64)

    # --------------------------
    # Helper: choose jobs sensibly
    # --------------------------
    def _n_jobs_for(n_perms: int) -> int:
        # don't spawn more workers than perms (wasteful)
        return max(1, min(int(n_cores), int(n_perms)))

    # --------------------------
    # Worker: run a chunk of permutations (for a dense chunk of genes)
    # Returns sums of scores and sums of squared scores for null mean/var
    # --------------------------
    def _perm_worker(X_chunk_centered, den_safe, n_perms, worker_seed):
        rng = np.random.default_rng(worker_seed)
        B = X_chunk_centered.shape[1]  # genes in this chunk
        sum_scores = np.zeros(B, dtype=np.float64)
        sum_sq = np.zeros(B, dtype=np.float64)
        L = lsa_graph  # closure

        for _ in range(n_perms):
            idx = rng.permutation(N)          # permute spot indices
            Xp = X_chunk_centered[idx, :]     # permute rows (spots)
            p_num = np.einsum("ij,ij->j", Xp, L @ Xp)  # numerator for each gene
            scores = p_num / den_safe
            sum_scores += scores
            sum_sq += scores * scores

        return sum_scores, sum_sq

    # --------------------------
    # Helper: run exactly n_perms_total permutations in parallel
    # --------------------------
    def _permute_chunk(X_chunk_centered, den_safe, n_perms_total, seed_base):
        n_jobs = _n_jobs_for(n_perms_total)
        base = n_perms_total // n_jobs
        rem = n_perms_total % n_jobs
        counts = [base + (1 if i < rem else 0) for i in range(n_jobs)]

        results = Parallel(n_jobs=n_jobs, prefer=prefer_backend)(
            delayed(_perm_worker)(
                X_chunk_centered,
                den_safe,
                counts[i],
                seed_base + i
            )
            for i in range(n_jobs) if counts[i] > 0
        )

        total_sum = np.sum([r[0] for r in results], axis=0)
        total_sum_sq = np.sum([r[1] for r in results], axis=0)
        return total_sum, total_sum_sq

    # --------------------------
    # PASS 1: observed scores + stage1 (or full) p-values for all genes
    # --------------------------
    if two_stage:
        n_pass1 = int(n_permutations_stage1)
        print(f"[Stage 1] {n_pass1} permutations for ALL genes")
    else:
        n_pass1 = int(n_permutations)
        print(f"[Single-stage] {n_pass1} permutations for ALL genes")

    print(f"Processing {G} genes in chunks of {chunk_size} for {N} spots/cells...")
    print(f"Backend: {prefer_backend} | dense dtype: {dtype}")

    for g_start in range(0, G, chunk_size):
        g_end = min(g_start + chunk_size, G)
        cols = slice(g_start, g_end)

        # Dense chunk: (N x B)
        X_chunk = X[:, cols].toarray().astype(dtype, copy=False)

        # Center each gene column by global mean
        X_chunk_centered = X_chunk - X_mean[cols].astype(dtype, copy=False)

        # Observed numerator/denominator
        spatially_smoothed = lsa_graph @ X_chunk_centered              # (N x B)
        num = np.einsum("ij,ij->j", X_chunk_centered, spatially_smoothed).astype(np.float64)
        den = np.sum(X_chunk_centered * X_chunk_centered, axis=0).astype(np.float64)
        den_safe = np.where(den == 0.0, 1e-12, den)

        obs_chunk = num / den_safe
        all_obs[g_start:g_end] = obs_chunk

        # Permutation moments for pass1
        seed_base = seed + (g_start // chunk_size) * 100000
        sum_scores, sum_sq = _permute_chunk(X_chunk_centered, den_safe, n_pass1, seed_base)

        mean_null = sum_scores / float(n_pass1)
        var_null = (sum_sq / float(n_pass1)) - mean_null * mean_null
        std_null = np.sqrt(np.maximum(var_null, 1e-12))

        z = (obs_chunk - mean_null) / std_null
        p_chunk = norm.sf(z)  # one-sided (positive spatial coherence)

        if two_stage:
            stage1_pvals[g_start:g_end] = p_chunk
        else:
            all_pvals[g_start:g_end] = p_chunk

        if g_start % 2000 == 0:
            print(f"Completed pass1: {g_end}/{G} genes...")

    # --------------------------
    # PASS 2: refine only top_k_refine genes (smallest stage1 p-values)
    # --------------------------
    if two_stage:
        # Start with stage1 p-values for all genes
        all_pvals[:] = stage1_pvals

        n_pass2 = int(n_permutations - n_permutations_stage1)
        if top_k_refine > 0 and n_pass2 > 0:
            # Choose genes to refine
            refine_idx = np.argsort(stage1_pvals)[:top_k_refine]
            refine_mask = np.zeros(G, dtype=bool)
            refine_mask[refine_idx] = True

            print(f"[Stage 2] +{n_pass2} perms for top_k_refine={top_k_refine} genes")

            for g_start in range(0, G, chunk_size):
                g_end = min(g_start + chunk_size, G)
                idx_chunk = np.arange(g_start, g_end)
                sel = refine_mask[idx_chunk]
                if not np.any(sel):
                    continue

                cols = idx_chunk[sel]  # actual gene indices we refine in this chunk

                # Dense chunk only for selected genes
                X_chunk = X[:, cols].toarray().astype(dtype, copy=False)
                X_chunk_centered = X_chunk - X_mean[cols].astype(dtype, copy=False)

                # Observed for selected genes
                spatially_smoothed = lsa_graph @ X_chunk_centered
                num = np.einsum("ij,ij->j", X_chunk_centered, spatially_smoothed).astype(np.float64)
                den = np.sum(X_chunk_centered * X_chunk_centered, axis=0).astype(np.float64)
                den_safe = np.where(den == 0.0, 1e-12, den)
                obs = num / den_safe

                # Recompute stage1 moments for this subset (cheap vs full G)
                seed_base1 = seed + (g_start // chunk_size) * 100000
                sum1, sq1 = _permute_chunk(X_chunk_centered, den_safe, n_permutations_stage1, seed_base1)

                # Stage2 additional moments
                seed_base2 = seed + 999999 + (g_start // chunk_size) * 100000
                sum2, sq2 = _permute_chunk(X_chunk_centered, den_safe, n_pass2, seed_base2)

                # Combine moments
                sum_scores = sum1 + sum2
                sum_sq = sq1 + sq2
                n_total = float(n_permutations_stage1 + n_pass2)

                mean_null = sum_scores / n_total
                var_null = (sum_sq / n_total) - mean_null * mean_null
                std_null = np.sqrt(np.maximum(var_null, 1e-12))

                z = (obs - mean_null) / std_null
                p_refined = norm.sf(z)

                all_pvals[cols] = p_refined

            print("[Stage 2] refinement complete.")
        else:
            print("[Stage 2] skipped.")

    # --------------------------
    # Final FDR across ALL genes
    # --------------------------
    _, fdr, _, _ = multipletests(all_pvals, alpha=alpha_svg, method="fdr_bh")

    # Store outputs
    adata.var["spartan_saq"] = all_obs
    adata.var["spartan_saq_pval"] = all_pvals
    adata.var["spartan_saq_fdr"] = fdr
    adata.var[key_added] = adata.var["spartan_saq_fdr"] < alpha_svg
    adata.var["spartan_saq_rank"] = (
        adata.var["spartan_saq"].rank(method="first", ascending=False).astype(int)
    )

    if copy:
        return adata
    return None
