KLDistance
===========

.. _KLDistance:

KL Distance Ranking
------------------

The KL (Kullback-Leibler) distance ranking function helps identify guides with significantly different spatial distributions compared to a reference distribution.

.. py:function:: rank_by_relative_entropy(adata, reference_guide='sum', control_guide='sgNon-targeting', result_field='KL distance', guide_list=None, n_top=50)

   Ranks guides by their KL divergence from a reference distribution.

   :param adata: AnnData object containing spatial transcriptomics data
   :param reference_guide: Reference distribution to compare against. Can be either 'sum' (sum of all guides) or 'ntc' (non-targeting control guide). Default is 'sum'.
   :param control_guide: Name of the non-targeting control guide in the dataset. Default is 'sgNon-targeting'.
   :param result_field: Name of the field to store results in adata.uns. Default is 'KL distance'.
   :param guide_list: Optional list of guides to analyze. If None, all guides will be analyzed. Default is None.
   :param n_top: Number of top guides to return in the results. Default is 50.
   :return: None. Results are stored in adata.uns[result_field] as a pandas DataFrame.

The function computes the KL divergence between each guide's spatial distribution and the reference distribution. A higher KL distance indicates a more distinct spatial pattern compared to the reference.

Example usage:

.. note:: 

    KL Distance Ranking is generally a method to model distribution of guides that have low spatial resolution or the spatial encoding is not the essential feature.
    As KL Distance Ranking dicards the spatial relationship between locations.

.. code-block:: 

    import B_cluster_independent
    B_cluster_independent.rank_by_relative_entropy(adata, reference_guide='sum')
    D_plot_ranking.plot_ranking(adata, 'KL distance')

The result is shown below.

.. image:: ../_images/KL_bar.png
   :align: center

We can check the distribution of the guides with high KL distance.

.. image:: ../_images/KL_hist.png
   :align: center

All KL distance results are stored in the :py:attr:`adata.uns` attribute named 'KL distance' by default.
