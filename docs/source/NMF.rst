NMF
====================

.. _NMF:

Gene program NMF Principle
---------------------------

**Principle:**  
Gene program NMF is a dimensionality reduction technique that factorizes a non-negative data matrix \( X \) into two lower-rank non-negative matrices \( W \) and \( H \), such that:

\[
X \approx W H
\]

where:
- \( X \) is the original data matrix (shape: \( n \times m \)), with all elements \( x_{ij} \geq 0 \),
- \( W \) (shape: \( n \times k \), \( k \ll m \)) is the basis matrix (or "bin"/"gene program" matrix),
- \( H \) (shape: \( k \times m \)) is the coefficient matrix (encodes the representation of original features in the lower-dimensional space).

**Optimization Problem:**  
Gene program NMF minimizes the difference between \( X \) and \( W H \) subject to non-negativity constraints:

\[
\min_{W, H} \ \|X - WH\|_F^2
\]
\[
\text{subject to: } W \geq 0, \ H \geq 0
\]

where \( \|\cdot\|_F \) indicates the Frobenius norm.

**Multiplicative Update Rules:**  
Gene program NMF uses a common algorithm for updating \( W \) and \( H \) is as follows (Lee and Seung, 2001):

\[
H_{kj} \leftarrow H_{kj} \cdot \frac{(W^T X)_{kj}}{(W^T W H)_{kj}}
\]
\[
W_{ik} \leftarrow W_{ik} \cdot \frac{(X H^T)_{ik}}{(W H H^T)_{ik}}
\]

**Interpretation:**  
- Each column of \( W \) defines a "gene program" or latent pattern.
- Each column of \( H \) describes the mixture coefficients (i.e., how much of each gene program is present) for the original samples/cells/spots.
- All matrices contain only non-negative values, enhancing interpretability for non-negative biological data (e.g., gene expression).

.. math::

   \mathbf{X} \approx \mathbf{W} \mathbf{H}
   
   \quad \text{with} \quad
   \mathbf{X} \in \mathbb{R}_{\geq 0}^{n \times m},
   \quad \mathbf{W} \in \mathbb{R}_{\geq 0}^{n \times k},
   \quad \mathbf{H} \in \mathbb{R}_{\geq 0}^{k \times m}

   \quad \Rightarrow \quad
   \underset{\mathbf{W},\,\mathbf{H}\,\geq\,0}{\mathrm{argmin}}
   \ \| \mathbf{X} - \mathbf{W} \mathbf{H} \|_F^2

Gene program NMF clustering using :py:func:`tardis_spac.utils.nmf_clustering()`

.. py:function:: tardis_spac.utils.nmf_clustering(adata, n_components=10, random_state=42, max_iter=1000, verbose=0, n_top_genes=2000)

   Perform NMF-based clustering on highly variable genes.

   :param adata: AnnData object. The raw expression matrix.
   :param n_components: Number of NMF components/clusters.
   :param random_state: Random seed for reproducibility.
   :param max_iter: Maximum number of NMF iterations.
   :param verbose: Verbosity level for fitting NMF.
   :param n_top_genes: Number of highly variable genes to use.
   :return: AnnData object containing NMF results in obsm['X_nmf'] and uns['X_nmf_components'].

.. admonition:: Usage

   With NMF clustering, highly variable genes that are 'dominantly' effecting the tumor landscape are identified.
   These genes are then used to cluster the spots into different gene programs.

Gene program NMF Consensus Clustering using :py:func:`tardis_spac.utils.nmf_consensus()`

.. py:function:: tardis_spac.utils.nmf_consensus(adata, min_clusters=4, max_clusters=10, n_resamples=100, resample_frac=0.8, random_state=42, n_cluster_genes=50)

   Perform consensus clustering using NMF results and compute cluster gene scores.

   :param adata: AnnData object where NMF results have already been computed.
   :param min_clusters: Minimum number of clusters to test.
   :param max_clusters: Maximum number of clusters to test.
   :param n_resamples: Number of resamplings for consensus clustering.
   :param resample_frac: Fraction of samples used in each resample.
   :param random_state: Random seed.
   :param n_cluster_genes: Number of top genes per cluster to use for scoring.
   :return: AnnData object with consensus clustering and normalized gene scores added to .obs.

.. admonition:: Usage

   With NMF consensus clustering, the gene programs are further refined by computing the consensus clustering of the NMF results.
   This is done by resampling the spots and computing the NMF results for each resample.
   The consensus clustering is then computed by averaging the NMF results across the resamples.
   The gene scores are then computed by averaging the NMF results across the resamples for each gene.
