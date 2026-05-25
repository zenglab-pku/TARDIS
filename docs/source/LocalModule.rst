LocalModule
===============

.. _LocalModule:

Guide Ranking by Cluster-Dependent Distribution
----------------------------------------------

This module provides two main statistical approaches for ranking guides in spatial transcriptomics data based on their distribution across clusters: **PERMANOVA-based ranking** and **Aitchison distance-based ranking**. Both methods compare each guide's distribution to that of a reference (usually a negative control guide, e.g., "sgNon-targeting") and support permutation testing for empirical significance assessment.

PERMANOVA-Based Clustering Guide Ranking
----------------------------------------

.. py:function:: permanova(gdata, cluster_field, result_field='permanova_f_value', reference_guide='sgNon-targeting', library_key=None, count_bins=10, n_permutations=1000, show_progress=True, copy=False)

   Ranks guides by using the PERMANOVA (Permutational Multivariate Analysis of Variance) test to compare their distributions across clusters to a specified reference guide. This method is adapted to spatial transcriptomics by comparing the distributions of cluster counts for each guide to the reference guide via kernel density estimation and the Bray-Curtis distance.

   **Mathematical Principle:**

   Let :math:`X_g = \{x_{gc}\}` be the vector of cell counts (or expression) for guide :math:`g` across each cluster :math:`c`. For guide :math:`g` and the reference guide, kernel density estimation is performed over the binned counts. The two group distributions are compared using Bray-Curtis dissimilarity. The PERMANOVA test pseudo F-statistic for the two groups is then computed as

   .. math::

      F = \frac{SS_B}{SS_W + \epsilon}

   where :math:`SS_B` is the between-group mean distance, :math:`SS_W` is the within-group mean distance, and :math:`\epsilon` is a small constant for numerical stability.

   Permutation testing is used to assess significance by randomly swapping cluster densities between the guide group and the reference across permutation runs.

   :param gdata: AnnData object with .obs (for cell metadata), .var (for guide info), and .X (expression/count matrix).
   :param cluster_field: Key in .obs indicating cluster labels.
   :param result_field: Field name in .var to write PERMANOVA F statistic (default: "permanova_f_value").
   :param reference_guide: The guide name used as the negative control/reference (default: "sgNon-targeting").
   :param library_key: If set, analysis will be performed separately per batch/sample group (default: None).
   :param count_bins: Number of bins for cluster-wise density estimation (default: 10).
   :param n_permutations: Number of permutations for empirical p-value estimation (default: 1000).
   :param show_progress: Whether to display progress during permutations (default: True).
   :param copy: If True, return a new modified AnnData object; else, modify in-place (default: False).

   :return: None if inplace; otherwise returns AnnData with results written to .var.

   Example:

   .. code-block::

      import tardis as td
      td.tardis_spac.stats.permanova(adata, cluster_field='leiden')
      adata.var.sort_values('permanova_f_value', ascending=False)

   After running, the field "permanova_f_value" in .var contains F statistics for each guide, and "permanova_f_value.p_value" contains empirical p-values from permutation tests.

   **Note:** The method compares the *shape* of cluster abundance distributions between guides and the reference, independently for each guide. Clustering assignment is required in advance.

Aitchison Distance-Based Guide Ranking
--------------------------------------

.. py:function:: aitchison_distance(gdata, cluster_field, result_field='aitchison_dist', reference_guide='sgNon-targeting', library_key=None, n_permutations=1000, show_progress=True, p_swap=0.1, copy=False)

   Computes the Aitchison distance between each guide and the reference guide based on their compositional (cluster-wise) abundances, ranking guides by the resulting distance. The method supports permutation testing by swapping cluster values between the guide and reference with probability `p_swap`.

   **Mathematical Principle:**

   Given the abundance composition of each guide :math:`g` over clusters :math:`c` (counts :math:`x_{gc}`), the composition vector is transformed as:

   .. math::

      v_{g} = \log(x_{gc} + 1) - \frac{1}{C} \sum_{c'} \log(x_{gc'} + 1)

   where C is the number of clusters.

   The Aitchison distance is then:

   .. math::

      d(g, r) = \sqrt{ \sum_c (v_{g,c} - v_{r,c})^2 }

   Permutation testing proceeds by swapping counts between guide and reference in each cluster with probability `p_swap` and recomputing Aitchison distances.

   :param gdata: AnnData object containing .obs, .var, .X.
   :param cluster_field: Key in .obs denoting cluster assignment for each cell.
   :param result_field: Output field name in .var for storing Aitchison distances (default: "aitchison_dist").
   :param reference_guide: Reference guide name (default: "sgNon-targeting").
   :param library_key: Key in .obs for performing per-sample (library) analysis (default: None).
   :param n_permutations: Number of permutations for p-value estimation (default: 1000).
   :param show_progress: Display progress bar for permutations (default: True).
   :param p_swap: Probability of swapping the cluster counts between the sample and reference per permutation (default: 0.1).
   :param copy: Return new AnnData if True; operate in-place if False (default: False).

   :return: None if inplace; otherwise modified AnnData.

   Example:

   .. code-block::

      import tardis as td
      td.tardis_spac.stats.aitchison_distance(adata, cluster_field='leiden')
      adata.var.sort_values('aitchison_dist', ascending=False)

   After running, "aitchison_dist" and "aitchison_dist.p_value" fields in .var will contain distance and corresponding empirical p-values.

   **Note:** This method treats guide cluster abundance as a *composition* and measures divergence from the reference using Aitchison geometry (Euclidean distance after log-ratio transform on composition). The permutation null ensures fair empirical significance control. Clustering assignment is required in advance.

Result Storage & Usage
----------------------

- All ranking and statistical results (scores and p-values) are written to fields in ``adata.var`` as specified by `result_field`.
- Use the result to filter, rank, or further visualize guides that drive distinct cluster distributions compared to controls.
- Both methods optionally allow per-sample or per-library group analysis via `library_key`, writing per-group results.
- Appropriate permutation-based p-values help control for statistical significance across the high-dimensional distribution space.
