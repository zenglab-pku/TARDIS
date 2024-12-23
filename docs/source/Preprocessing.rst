Preprocessing
=====

.. _Preprocessing:

.. note::

   This here is an API version of all the functions applicable in **SPOT**,
   For more detailed and thorough reference please visit our :ref:`tutorial` site.

Preprocessing
----------------

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
   
