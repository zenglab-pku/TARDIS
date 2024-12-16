Preprocessing
=====

.. _Preprocessing:

.. note::

   This here is an API version of all the functions applicable in **SPOT**,
   For more detailed and thorough reference please visit our :ref:`tutorial` site.

Preprocessing
----------------

Spatial perturbation can be highly arbitrary if we cannot perform valid
preprocessing and filtering of low quality guides and bins. Refer to *Paper name*
for our in house filtering method.

**SPOT** performs filtering with validation panels with the following methods.

First import **SPOT** from library.

.. code-block::

   import spot as sp
   sp.set_random_seed(42)

In this tutorial we use our in house spatial transcriptomics data.
This data incorporates a library of **68 guides**, and is sequenced on **BGI Stereo-seq** platform.

.. code-block::

   # perform quality check from BGI stereo-seq GEM output
   sp.preprocessing.filter_qc_bins('A04091E1.gem')

*output:*

.. image:: ../_images/qc_guide_bins.png
   :align: center

Filtering using :py:func:`preprocessing.filter_guide_reads()`

.. py:function:: filter_guide_reads(gem_path, guide_prefix=None, output_path=None, binarilize=False, assign_pattern='max', filter_threshold=None)

   Filter and process guide reads from a GEM file.

   :param gem_path: Path to the input GEM file
   :param guide_prefix: Optional prefix to filter guide names. Only guides starting with this prefix will be kept
   :param output_path: Optional path to save filtered results. If None, returns the filtered DataFrame
   :param binarilize: Whether to set all bin counts to 1. Recommended for high library size and resolution
   :param assign_pattern: How to handle bins with multiple guides. Can be 'max' (keep guide with max count), 'drop' (remove multi-guide bins), or 'all' (keep all guides)
   :param filter_threshold: Minimum guide count threshold. Bins with guides below this count will be filtered out
   :return: Filtered pandas DataFrame if output_path is None, otherwise None

The function processes a GEM file containing guide reads and performs filtering based on the specified parameters:

1. Reads the GEM file and optionally filters for guides with a specific prefix
2. Removes bins with guide counts below the threshold if specified  
3. Handles bins with multiple guides according to the assign_pattern:

   - 'max': Keeps only the guide with highest count in each bin
   - 'drop': Removes all bins that have multiple guides
   - 'all': Keeps all guides in multi-guide bins

4. Optionally binarizes the counts (sets all to 1)
5. Returns filtered DataFrame or saves to file

Example usage:

.. code-block::

   filtered_data = sp.filter_guide_reads('A04091E1.gem', output_path='A04091E1_filtered.gem')

After filtering, we can perform quality control on the filtered data.

.. code-block::

   sp.preprocessing.filter_qc_bins('A04091E1_filtered.gem')

   plt.figure(figsize=(8, 6))
   scatter = plt.scatter(x=gdata.obsm['spatial'][:, 0], y=gdata.obsm['spatial'][:, 1],
                        s=gdata.obs['n_genes_by_counts'], alpha=0.5, c=gdata.obs['total_counts'], cmap='viridis')
   sns.despine()
   plt.colorbar(scatter, label='Total counts')
   plt.title('Guide reads')

   plt.show()

.. image:: ../_images/guide_reads.png
   :align: center
