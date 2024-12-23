NMF
====================

.. _NMF:

.. py:function:: A_preprocessing.nmf_clustering(adata, n_components=10, random_state=42, max_iter=1000, verbose=0, n_top_genes=2000)

   Performs NMF clustering on AnnData object.

   :param adata: AnnData object containing gene expression data
   :param n_components: Number of NMF components to extract
   :param random_state: Random seed for reproducibility
   :param max_iter: Maximum number of iterations for NMF
   :param verbose: Verbosity level
   :param n_top_genes: Number of highly variable genes to use
   :return: AnnData object with NMF results in .obsm['X_nmf'] and .uns['X_nmf_components']

.. py:function:: A_preprocessing.nmf_consensus(adata, min_clusters=4, max_clusters=10, n_resamples=100, resample_frac=0.8, random_state=42, n_cluster_genes=50)

   Performs consensus clustering on NMF results.

   :param adata: AnnData object with NMF results
   :param min_clusters: Minimum number of clusters to try
   :param max_clusters: Maximum number of clusters to try  
   :param n_resamples: Number of resampling iterations
   :param resample_frac: Fraction of samples to use in each resample
   :param random_state: Random seed for reproducibility
   :param n_cluster_genes: Number of top genes to use per cluster
   :return: AnnData object with cluster scores in .obs
