IO
========

.. _IO:

This section describes the input of **TARDIS**.

Loading Stereo-seq Data
----------------------

.. py:function:: tardis.io.load_stereo

    .. code-block::

        import tardis as td
        bin_df = td.io.load_bin('/path/to/your/guide.gem',
                                bin_size=100,
                                library_id='spatial')

    :param path: str
        The path to the Stereo-seq data.
    :param library_id: str
        The library id to save in obsm['spatial'].
    :param bin_size: int
        The size of the bin. (Default: 100)

Loading Visium HD Data
----------------------

.. py:function:: tardis.io.load_visium

    .. code-block::

        import tardis as td
        bin_df = td.io.load_visium('/path/to/your/visium_hd_data',
                                   library_id='spatial',
                                   use_array=False)

    :param path: str
        The path to the Visium HD data.
    :param library_id: str
        The library id to save in obsm['spatial'].
    :param use_array: bool
        Whether to use the array row and col as coordinates.
