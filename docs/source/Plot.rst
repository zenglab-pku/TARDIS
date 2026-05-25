Plotting Functions
==================

.. _Plot:

TARDIS Plotting API
-------------------

TARDIS provides convenient plotting functions for visualizing guide-level statistics, top-ranked guides, and their spatial distribution.

These are useful for visually inspecting the results of Aitchison distance or PERMANOVA-based guide ranking, and for spatial mapping of perturbation effects.

Plot KDE of Ranking Results
---------------------------

.. py:function:: tardis_spac.utils.plot_top_kde(gdata, result_field, show_sgnt=True, top_n=2, violin=True, sgnt_label="sgNon-targeting", figsize=(5,2.7), min_count=0)

   Plot the Kernel Density Estimate (KDE) distribution of a result field (e.g. Aitchison distance or PERMANOVA p-value) across all guides. Optionally highlight the sgNon-targeting guide and top N guides.

   :param gdata: AnnData object containing guide-level statistics in ``.var`` (e.g. result of ranking).
   :param result_field: Name of the column in ``gdata.var`` holding scores to plot.
   :param show_sgnt: Whether to highlight the sgNon-targeting guide. Default ``True``.
   :param top_n: Number of top guides to label and highlight. Default ``2``.
   :param violin: Whether to plot 'violin' vertical lines for each data point. Default ``True``.
   :param sgnt_label: Index (name) in ``gdata.var`` for the sgNon-targeting guide. Default ``"sgNon-targeting"``.
   :param figsize: Tuple. Size of plot in inches. Default ``(5,2.7)``.
   :param min_count: Minimum value in ``TotalCount`` for a guide to be included. Default ``0``.
   :returns: Matplotlib figure and axes.

   Use this function after running a ranking procedure on your data, such as :py:func:`rank_by_aitchison_distance`, to visualize the distribution of guide statistics.

   Example usage:

   .. code-block:: python

      import tardis_spac.utils as tu
      tu.plot_top_kde(adata, result_field="Aitchison distance", top_n=3)

   The sgNon-targeting control and top-ranking guides will be annotated on the KDE plot for easy interpretation.

Spatial Plotting of Guides
--------------------------

.. py:function:: tardis_spac.utils.plot_spatial_guides(gdata, scale_factor, image=None, figsize=(10,10), palette='tab20b', s=3, edgecolor='none', legend=False, alpha=1.0, ax=None)

   Visualize the spatial distribution of guides on a tissue section or reference image.

   :param gdata: AnnData object with guide count data; must have coordinates in ``.obsm['spatial']`` and guides as ``.var_names``.
   :param scale_factor: Float. Scaling to match spot coordinates to image pixels.
   :param image: Optional. Background image (`numpy.ndarray` or PIL Image) for the tissue section.
   :param figsize: Tuple giving figure size. Default ``(10,10)``.
   :param palette: Color palette for guides. Default ``'tab20b'``.
   :param s: Dot size. Default ``3``.
   :param edgecolor: Dot edge color. Default ``'none'``.
   :param legend: Whether to include the legend. Default ``False``.
   :param alpha: Transparency level for dots. Default ``1.0``.
   :param ax: Matplotlib axis to plot on (optional).
   :returns: Matplotlib axes.

   Use this function to visualize where each guide is spatially enriched or distributed across the tissue. Dots are colored by the most abundant guide for each spatial bin.

   Example usage:

   .. code-block:: python

      import tardis_spac.utils as tu
      tu.plot_spatial_guides(
          gdata, 
          scale_factor=0.5, 
          image=tissue_img, 
          palette='tab20', 
          s=4
      )

   This will produce a plot overlaying guide locations on the provided image and coloring by assigned guide identity.

Tips and Interpretation
-----------------------

- Use :py:func:`plot_top_kde` after guide-level ranking (Aitchison or PERMANOVA) to identify top outlier guides compared to the control.
- :py:func:`plot_spatial_guides` is ideal for checking the localization and targeting effect of guides in spatial-omics data.
- Adjust ``scale_factor`` so that spot coordinates match the scale of your background tissue image if used.
- Plots are returned as Matplotlib objects, so you can further customize or save them as needed.

.. note::

   Requires ``matplotlib``, ``numpy``, ``seaborn``, and ``statsmodels`` (for KDE). For further customization, modify the plotting code or provide additional Matplotlib axes as needed.

See Also
--------

- :doc:`Tutorial`
