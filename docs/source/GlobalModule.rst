GlobalModule
=================

.. _GlobalModule:

Kernel Density Estimation and Wasserstein Distance
------------------------------------------------------

Overview
--------

The functions in this module provide statistical ranking of guides in spatial transcriptomics data by quantifying the difference between the spatial distributions induced by each guide and a control, using nonparametric statistics and information theory. These metrics are essential for identifying guides that significantly alter spatial gene expression patterns, pointing to candidate genes that may create or disrupt spatial cellular "niches" — microenvironments with specific cellular compositions or molecular signaling cues. 

Biological Rationale
--------------------
In complex tissues, cellular functions and states are influenced not only by gene expression but also by spatial context—the "niche" a cell is in. By ranking perturbations (guides) based on how much they reshape spatial patterns (either absolutely using the KL divergence, or with explicit spatial information via the Wasserstein distance), the methods here help identify genes crucial for maintaining or disrupting these spatial niches.

Mathematical Principles
-----------------------

**Kernel Density Estimation (KDE):**
  - KDE provides a smooth estimate of the spatial probability density function of the transcriptomic signal for each guide. 
  - Mathematically, given spatial points \( x_1, x_2, ..., x_n \), the KDE at a location \( x \) is:
    \[
    \hat{f}(x) = \frac{1}{nh} \sum_{i=1}^n K\left( \frac{x - x_i}{h} \right)
    \]
    where \( K \) is a kernel function (e.g., Gaussian), and \( h \) is the bandwidth.
    
**Wasserstein (Earth Mover's) Distance:**
  - The Wasserstein distance quantifies the minimal "cost" to transform one spatial distribution into another, reflecting both spatial and distributional differences.
  - In 1D, for probability distributions \( u \) and \( v \):
    \[
    W(u, v) = \inf_{\gamma \in \Gamma(u, v)} \int |x-y| d\gamma(x, y)
    \]
    where \( \Gamma \) denotes all joint distributions with marginals \( u \) and \( v \).
  - It is sensitive to both the magnitude and the locations of expression, making it ideal for spatial data.

Permutation Significance:
  - For both statistics, statistical significance (empirical p-value) is assessed by random permutations of guide labels; the p-value is the fraction of permutations yielding a more extreme value than observed.

**KL (Kullback-Leibler) Divergence:**
  - KL divergence measures the relative entropy between two discrete distributions. For probability distributions \( P \) and \( Q \):
    \[
    D_{KL}(P \| Q) = \sum_i P(i) \log \frac{P(i)}{Q(i)}
    \]
  - In this context, it quantifies how much the overall expression pattern of a guide differs from a reference (e.g., a non-targeting guide).

Function Documentation
----------------------

.. py:function:: wasserstein_distance(adata, control_guide='sgNon-targeting', guide_list=None, n_permutation=50, n_process=8, return_fig=False, return_dataframe=True, sort_by_replicate='_')

   Ranks guide perturbations by how much they change the spatial distribution of expression, as measured by kernel density estimation and Wasserstein distance from a control guide.

   :param adata: AnnData object with spatial transcriptomics data, `.obsm[spatial]` specifying coordinates.
   :param control_guide: Name of the reference or negative control guide. (default: "sgNon-targeting")
   :param guide_list: List of guides to analyze (default: None, uses all guides except control_guide)
   :param n_permutation: Number of permutations for empirical p-value estimation. (default: 50)
   :param n_process: Number of parallel processes for permutation tests. (default: 8)
   :param return_fig: If True, returns Matplotlib figure. (default: False)
   :param return_dataframe: If True, returns results as DataFrame. (default: True)
   :param sort_by_replicate: Delimiter for identifying replicates (default: "_")

   :return: Depending on arguments, a result DataFrame, a figure, or both.

   **Method**
      1. For each guide, estimate its spatial distribution (KDE).
      2. Compute the Wasserstein distance from this guide to the control.
      3. Estimate empirical p-value for observed distance by comparing to a null distribution from permutations.
      4. Higher Wasserstein distances (with significant p-value) suggest a guide creates a new spatial niche or disrupts existing ones.

   Example:

   .. code-block::

      import tardis as td
      results = td.stats.wasserstein_distance(adata, control_guide='sgNon-targeting')

   **Interpretation:** Guides ranked at the top most strongly perturb spatial structure. A strong, significant Wasserstein distance means the guide changes the geography of gene expression, pointing to niche-altering genes.

KL Divergence-Based Ranking
---------------------------

.. py:function:: kl_divergence(adata, reference_guide='sum', control_guide='sgNon-targeting', result_field='KL distance', guide_list=None, n_top=50)

   Ranks guides by the Kullback-Leibler (KL) divergence of their expression distributions relative to a specified reference.

   :param adata: AnnData object.
   :param reference_guide: Reference distribution to use ("sum" for aggregate/background spatial profile, or a specific control guide such as "sgNon-targeting"). (default: "sum")
   :param control_guide: Name of the non-targeting or negative control guide. (default: "sgNon-targeting")
   :param result_field: Field name in `adata.uns` for storing results. (default: "KL distance")
   :param guide_list: Guides to analyze (default: None, uses all guides)
   :param n_top: Number of top guides to include in ranking. (default: 50)

   :return: None (results are stored in `adata.uns[result_field]` as a pandas DataFrame).

   **Method**
      1. Normalize each guide's total expression as a probability distribution across cells or bins.
      2. Compute the KL divergence from each guide's distribution to the reference.
      3. High values suggest the guide induces a distinct expression pattern (but not necessarily a spatially localized "niche"—see Note below).

   Example:

   .. code-block:: 

      import tardis as td
      td.stats.kl_divergence(adata, reference_guide='sum')
      # To visualize: td.plot_ranking.plot_ranking(adata, 'KL distance')

   .. note:: 

      KL divergence-based ranking is best suited for cases where spatial encoding is not central (e.g., sparse or low-resolution spatial data, or when distinguishing differences in global expression profile rather than explicit spatial localization). KL divergence considers the distribution across locations but discards their physical spatial relationships.

**Storage of Results**

- All Wasserstein distance and KL divergence results, including statistical significance and rankings, are stored as columns in ``adata.var`` (and, for summary tables, in ``adata.uns``), with keys such as ``w_dist``, ``w_dist.p_value``, ``KL distance``, etc.

**References**
  - [1] Rubner, Y., Tomasi, C., & Guibas, L. J. (2000). The earth mover's distance as a metric for image retrieval. International Journal of Computer Vision.
  - [2] Kullback, S., & Leibler, R. A. (1951). On information and sufficiency. Annals of Mathematical Statistics.
  - [3] Schiebinger, G., et al. (2019). Optimal-transport analysis of single-cell gene expression identifies developmental trajectories in reprogramming. Cell.

**Biological Interpretation**
  - Guides with high and significant distance values may define, induce, or disrupt unique cellular neighborhoods ("niches") within the tissue, shedding light on the molecular mechanisms underlying spatial organization.
