NMF
====================

.. _NMF:

NMF Clustering using :py:func:`tardis.preprocessing.nmf_clustering()`

.. py:function:: tardis.preprocessing.nmf_clustering(adata, n_components=10, random_state=42, max_iter=1000, verbose=0, n_top_genes=2000)

   Perform NMF-based clustering on highly variable genes.

   :param adata: AnnData object. The raw expression matrix.
   :param n_components: Number of NMF components/clusters.
   :param random_state: Random seed for reproducibility.
   :param max_iter: Maximum number of NMF iterations.
   :param verbose: Verbosity level for fitting NMF.
   :param n_top_genes: Number of highly variable genes to use.
   :return: AnnData object containing NMF results in obsm['X_nmf'] and uns['X_nmf_components'].

NMF Consensus Clustering using :py:func:`tardis.preprocessing.nmf_consensus()`

.. py:function:: tardis.preprocessing.nmf_consensus(adata, min_clusters=4, max_clusters=10, n_resamples=100, resample_frac=0.8, random_state=42, n_cluster_genes=50)

   Perform consensus clustering using NMF results and compute cluster gene scores.

   :param adata: AnnData object where NMF results have already been computed.
   :param min_clusters: Minimum number of clusters to test.
   :param max_clusters: Maximum number of clusters to test.
   :param n_resamples: Number of resamplings for consensus clustering.
   :param resample_frac: Fraction of samples used in each resample.
   :param random_state: Random seed.
   :param n_cluster_genes: Number of top genes per cluster to use for scoring.
   :return: AnnData object with consensus clustering and normalized gene scores added to .obs.
