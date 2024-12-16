API
===

.. _API:

API functions
------------------

.. py:function:: load_bin(gem_file, bin_size, library_id, image_file=None)

   Read BGI data and image file, and return an AnnData object.

   :param gem_file: The path of the BGI data file
   :param bin_size: The size of the bin
   :param library_id: The library id
   :param image_file: The path of the image file (optional)
   :return: AnnData object with the following keys:

      - :attr:`anndata.AnnData.obsm` ``['spatial']`` - spatial spot coordinates
      - :attr:`anndata.AnnData.uns` ``['spatial']['{library_id}']['images']`` - *hires* images
      - :attr:`anndata.AnnData.uns` ``['spatial']['{library_id}']['scalefactors']`` - scale factors for the spots


.. autosummary::
   :toctree: generated

   lumache
