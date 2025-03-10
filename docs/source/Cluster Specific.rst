Cluster Specific
=================

.. _ClusterSpecific:


To account for spatially clustered perturbations forming distinct clones, as observable in Lung metastatis data.
Here the bottom-up mode of TARDIS is recommended.
This method identifies perturbation-induced clonal formations while preserving their spatial specificity.
This is particularly relevant when perturbations drive localized tumor microenvironment changes that may be diluted if grouped into predefined micro-niches.

Clonal Defenition
-----------------

.. py:function:: clonal_definition(adata, cluster_field, control_guide='sgNon-targeting', n_permutation=999, guide_list=None, result_field='PERMANOVA p-value')

   Ranks guides by PERMANOVA test comparing their cluster distributions to a control guide.

   :param adata: AnnData object containing spatial transcriptomics data
   :param cluster_field: Name of the clustering field in adata.obs to use for grouping
   :param control_guide: Name of the control guide to compare against. Default is 'sgNon-targeting'.
   :param n_permutation: Number of permutations for PERMANOVA test. Default is 999.
   :param guide_list: Optional list of guides to analyze. If None, all guides will be analyzed. Default is None.
   :param result_field: Name of the field to store results in adata.uns. Default is 'PERMANOVA p-value'.
   :return: None. Results are stored in adata.uns[result_field] as a pandas DataFrame.

The function performs PERMANOVA tests comparing each guide's distribution across clusters to the control guide's distribution. A lower p-value indicates a more significant difference from the control guide pattern.

Example usage:

.. code-block::

    import tardis as td
    td.cluster_specific.clonal_definition(adata, cluster_field='leiden')
    td.cluster_specific.plot_ranking(adata, 'PERMANOVA p-value')

.. note::

    Clonal definition requires clustering information. Make sure to perform clustering on your data first.
    The test evaluates whether guides show significantly different patterns across clusters compared to the control using permutation-based statistical testing.
