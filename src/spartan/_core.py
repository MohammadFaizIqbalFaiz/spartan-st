from __future__ import annotations

from typing import Optional, Tuple, Union
import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix

import scanpy as sc
import squidpy as sq
from anndata import AnnData

import igraph as ig
import leidenalg
from joblib import Parallel, delayed

from scipy.stats import norm
from statsmodels.stats.multitest import multipletests


def normalize(mat): 
   row_means = np.array(mat.sum(axis=1)).flatten()#Calculates mean of each row
   row_means[row_means==0]=1#Avoids zero means
   inv_rm = sp.diags(1.0/row_means)#Creates sparse diagonal matrix where each diagonal element is equal to the row mean
   mat = inv_rm.dot(mat)#Calculates dot product of Spatial weight matrix and the diagonal matrix

   return mat

#robust vector norm

def vec_norm(v1,v2,flag:str):
    

    if(v1.shape!=v2.shape):
        print("vectors don't align")
    else:
        vt = v1 - v2
        #vt = np.maximum(vt,0.0)
   
    if flag == 'N':
        return np.linalg.norm(vt)
    elif flag == 'S':
        return np.linalg.norm(vt,axis=1)

def build_spatial_graph(adata:AnnData,crd_type:str,neighs:int,rings:int,neighborhood:str):
    
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
    
    #Perform PCA

    sc.pp.scale(adata)#Centering of the data
    sc.pp.pca(adata, n_comps=comps, svd_solver='arpack',random_state = seed)#Performs PCA using only highly variable genes

    PCA = adata.obsm['X_pca'][:, :comps_extract]#Extracts PCA values

    return adata,PCA

def build_lsa_graph(adata:AnnData,expr_values,spatial_matrix):

   
   G = spatial_matrix.tocoo()#Converts it into coco format

    #Calculates local mean gene expression of each neighborhood
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
    
    
    
   row = G.row
   col = G.col
   data = G.data

   xi = expr_values[row]
   xj= expr_values[col]

   local_mean_i = local_means[row]
   local_var_i = local_vars[row]

   diffI = vec_norm(xi,local_mean_i,flag='S')
   diffI = np.nan_to_num(diffI, nan=0.0)
   diffJ = vec_norm(xj,local_mean_i,flag='S')
   diffJ = np.nan_to_num(diffJ, nan=0.0)

   attribute_deviation = diffI * diffJ
   attribute_deviation = attribute_deviation / (local_var_i + np.finfo(float).eps)
   adjusted_weight = data * attribute_deviation
   G = coo_matrix((adjusted_weight, (row, col)), shape=G.shape)

   G = G.tocsr()

   return G

def build_gene_expression_connectivity_graph(adata:AnnData,crd_type:str,neighs:int):

    
    #Create neighborhood graphs based on gene expression similarity in PCA coordinate space
    sq.gr.spatial_neighbors(adata, coord_type=crd_type, n_neighs=neighs, key_added = 'gene_graph')
    gene_graph = adata.obsp['gene_graph_connectivities']

    return gene_graph

def to_igraph(adjacency):
    sources, targets = adjacency.nonzero()
    weights = adjacency[sources, targets].A1  # Convert sparse matrix to flat array
    g = ig.Graph(directed=False)
    g.add_vertices(adjacency.shape[0])
    g.add_edges(zip(sources, targets))
    g.es['weight'] = weights
    return g

def perform_leiden_clustering(adata:AnnData,lsa_graph:csr_matrix,spatial_graph:csr_matrix,gene_graph:csr_matrix,
                              alpha:float,beta1:float,beta2:float,resolution:float,seed:int,key_added:str):

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

def spartan_spatial_domains(adata:AnnData,spatial_coord:str="grid",spatial_neighs:Optional[int]=6,spatial_rings:Optional[int]=2,
                            spatial_neighborhood:str="knn",
                            total_pca_comps:int=50,pca_comps_extract:int=30,gene_coord:str="generic",gene_neighs:int=15,
                            alpha:float=0.80,beta1:float=0.1,beta2:float=0.4,resolution:float=1.0,seed:int=1,key_added="spartan_domains",
                            copy: bool=False):
    if copy:
        adata=adata.copy()
    
    adata,spatial_graph,spatial_W = build_spatial_graph(adata,crd_type=spatial_coord,neighs=spatial_neighs,rings=spatial_rings,
                                                       neighborhood=spatial_neighborhood)
    adata.obsp['spartan_spatial_graph'] = spatial_graph
    adata.obsp['spartan_spatial_weights'] = spatial_W
    
    adata,PCA = perform_pca(adata,comps=total_pca_comps,comps_extract=pca_comps_extract,seed=seed)

    expression_values = PCA.toarray() if hasattr(PCA, 'toarray') else PCA#Creates a dense array

    lsa_graph = build_lsa_graph(adata,expr_values=expression_values,spatial_matrix=spatial_W)

    adata.obsp['spartan_lsa_graph'] = lsa_graph

    #construct gene expression connectivity graph
    spac = adata.obsm['spatial']
    adata.obsm['spatial'] = PCA #Uses PCA coordinates to create gene expression similarity based neighborhood graphs

    gene_graph = build_gene_expression_connectivity_graph(adata,crd_type=gene_coord,neighs=gene_neighs)

    adata.obsp['spartan_gene_graph']=gene_graph
    adata.obsm["spatial"]=spac
    adata = perform_leiden_clustering(adata,lsa_graph,spatial_graph,gene_graph,alpha,beta1,beta2,resolution,seed,key_added)

    if copy:
        return adata
    return None




def spartan_svg(adata:AnnData, lsa_graph:csr_matrix, layer: str | None = "log1pX",
                n_permutations:int=1000, n_cores:int=8, use_X_if_missing: bool = True,
                                 alpha_svg:float =0.05, chunk_size:int=200, seed:int=1,key_added="spartan_svg",copy: bool=False):
    """
    Optimized for 80k+ cells and 19k+ genes.
    Uses chunking to stay within 128GB RAM.
    """
    if copy:
        adata=adata.copy()
    if layer is not None and layer in adata.layers:
        X = adata.layers[layer]
    elif use_X_if_missing:
        X = adata.X
    else:
        raise KeyError(
            f"Layer {layer!r} not found in adata.layers "
            "and use_X_if_missing=False."
        )    
    # 1. Ensure X is sparse and lsa_graph is CSR
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    if not sp.issparse(lsa_graph):
        lsa_graph = sp.csr_matrix(lsa_graph)
    
    N, G = X.shape
    all_obs = np.zeros(G)
    all_pvals = np.ones(G)

    # Pre-calculate global stats to avoid doing it per permutation
    # Mean of sparse matrix across cells (columns)
    X_mean = np.array(X.mean(axis=0)).flatten()
    
    print(f"Processing {G} genes in chunks of {chunk_size} for {N} cells...")

    # Process in chunks of genes
    for g_start in range(0, G, chunk_size):
        g_end = min(g_start + chunk_size, G)
        
        # Pull out a dense chunk: (Cells x Chunk_Size)
        # This is where we control the RAM
        X_chunk = X[:, g_start:g_end].toarray()
        X_chunk_centered = X_chunk - X_mean[g_start:g_end]
        
        # Observed Score for Chunk
        # (N x N) sparse @ (N x Chunk) dense
        spatially_smoothed = lsa_graph @ X_chunk_centered
        num = np.einsum("ij,ij->j", X_chunk_centered, spatially_smoothed)
        den = np.sum(X_chunk_centered**2, axis=0)
        den_safe = np.where(den == 0, 1.0, den)
        obs_chunk = num / den_safe
        all_obs[g_start:g_end] = obs_chunk

        # Parallel Permutations for this chunk
        def run_perms(n_p, s):
            rng = np.random.default_rng(s)
            s_p = np.zeros(g_end - g_start)
            sq_p = np.zeros(g_end - g_start)
            for _ in range(n_p):
                idx = rng.permutation(N)
                Xp = X_chunk_centered[idx, :]
                p_num = np.einsum("ij,ij->j", Xp, lsa_graph @ Xp)
                scores = p_num / den_safe
                s_p += scores
                sq_p += scores**2
            return s_p, sq_p

        perms_per_core = n_permutations // n_cores
        results = Parallel(n_jobs=n_cores)(
            delayed(run_perms)(perms_per_core, seed + i) for i in range(n_cores)
        )

        # Aggregate Chunk Results
        total_sum = np.sum([r[0] for r in results], axis=0)
        total_sum_sq = np.sum([r[1] for r in results], axis=0)
        
        m_null = total_sum / n_permutations
        v_null = (total_sum_sq / n_permutations) - m_null**2
        s_null = np.sqrt(np.maximum(v_null, 1e-12))
        
        z = (obs_chunk - m_null) / s_null
        all_pvals[g_start:g_end] = norm.sf(z)
        
        if g_start % 1000 == 0:
            print(f"Completed {g_end}/{G} genes...")

    # Final FDR
    _, fdr, _, _ = multipletests(all_pvals, alpha=alpha_svg, method='fdr_bh')
    #significant = fdr < alpha_svg
    adata.var['spartan_saq'] = all_obs
    adata.var['spartan_saq_pval'] = all_pvals
    adata.var['spartan_saq_fdr'] = fdr
    adata.var[key_added] = adata.var['spartan_saq_fdr'] < alpha_svg
    adata.var['spartan_saq_rank'] = (
        adata.var['spartan_saq']
        .rank(method="first",ascending=False)
        .astype(int)
    )
    if copy:
        return adata
    return None
