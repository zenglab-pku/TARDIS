Preprocessing
=====

.. _Preprocessing:

.. note::

   This here is an API version of all the functions applicable in **TARDIS**,
   For more detailed and thorough reference please visit our :ref:`tutorial` site.

Visium HD
---------

With Visium HD data, you should run `Spaceranger <https://www.10xgenomics.com/support/software/space-ranger/latest>`_ with your guide reference and probe set reference.

Try the following tutorial for **Visium HD** customized reference building:
- `Customized guide reference building <https://www.10xgenomics.com/support/software/space-ranger/latest/advanced/custom-references>`_

Successfully setup, you should be able to find guide reads in your *filtered_bc_matrix.mtx* file in Spaceranger output directory.

Read Spaceranger output
~~~~~~~~~~~~~~~~~~~~~~~

You can simply apply `scanpy.read_10x_mtx <https://scanpy.readthedocs.io/en/stable/api/scanpy.read_10x_mtx.html>`_ to read the guide reads in your *filtered_bc_matrix.mtx* file.

Remember to set genome only to **FALSE** for proper guide read assignment.

Stereoseq
---------

With BGI stereoseq data, you should follow the instructions on `SAW <https://github.com/STOmics/SAW>`_ to process the data.

Use **SPAC-seq** script *SeekSeq* to process your enriched guide FASTQ files.

.. note::

   **SeekSeq** is a tool for processing stereoseq data.
   It is not a part of **TARDIS** toolkit.
   refer to **SPAC-seq** official code repository `SPAC-seq <https://github.com/zenglab-pku/SPAC-seq/blob/main/fig_3_seekseq.py>`_ for more details.

   In **SeekSeq**, we filter the guide reads by the following criteria:
   - R1 containing constant region "TTGTCTTCCTAAGAC" with a max error rate of 1 base mismatch.
   - R2 containing scaffold region "GTTTTAGA" with a max error rate of 1 base mismatch.

   See our paper for more details: `Uncovering Spatially Resolved Functional Genomics with CRISPR Screen Sequencing <https://www.cell.com/cell/fulltext/S0092-8674(26)00516-7>`_..

Successfully processed, you should be able to find guide reads in your *guide.gem* file in SAW output directory.

Read SAW output
~~~~~~~~~~~~~~~

Read in SAW output from BGI stereoseq data with :py:func:`tardis_spac.utils.load_bin()`

.. py:function:: tardis_spac.utils.load_bin(gem_file, bin_size)

   Read in SAW output from BGI stereoseq data with :py:func:`tardis_spac.utils.load_bin()`

   :param gem_file: Path to the input SAW GEM file.
   :param bin_size: The size of the bin.
   :return: AnnData object with the guide reads mapped to the bin.

.. note::
   For Stereoseq data only.

.. admonition:: Function Explanation

   The :py:func:`tardis_spac.utils.load_bin` function reads a SAW GEM file from BGI Stereoseq data and aggregates UMI (Unique Molecular Identifier) counts into spatial bins of a specified size. It generates an AnnData object whose spots represent these bins, making downstream analysis (such as integration with gene expression/transcriptome data) more convenient.

   **Choosing the right ``bin_size`` is important:**  
   - A small ``bin_size`` provides higher spatial resolution but may result in many bins with few reads. (20-50um)
   - A large ``bin_size`` increases read counts per bin but reduces spatial resolution. (50-100um)
   - For integrative analysis (e.g., with spatial transcriptome data), select a ``bin_size`` consistent with the transcriptome binning.

   **Tip**: You can read and bin the guide GEM file and transcriptome GEM file together using the same ``bin_size`` for compatible spatial mapping.

   **Tip**: For BGI stereoseq data, a bin is approximately 500nm in size.

.. code-block:: python

   from tardis_spac.utils import load_bin
   guide_adata = load_bin("guide.gem", bin_size=50)
   gene_adata = load_bin("transcriptome.gem", bin_size=50)
   # Now guide_adata and gene_adata are spatially comparable

TARDIS preprocessing functions
------------------------------

Quality Control of Guide Bins using :py:func:`tardis_spac.utils.qc_guide_bins()`

.. note::
   For Stereoseq data only.
   For Visium HD data, you can apply `scanpy.pp.calculate_qc_metrics <https://scanpy.readthedocs.io/en/stable/api/scanpy.pp.calculate_qc_metrics.html>`_ to calculate the quality metrics of the guide bins.

.. py:function:: tardis_spac.utils.qc_guide_bins(gem_path, guide_prefix=None, fig_path=None)

   Summarize and visualize the quality of guide bins in the provided GEM file.

   :param gem_path: Path to the input GEM file.
   :param guide_prefix: Optional string. If provided, filter guide IDs starting with this prefix.
   :param fig_path: Optional string. Path to save the output plot image. If None, the plot will be displayed instead.
   :return: None. Shows or saves a figure displaying summary statistics of the guide bins.

Remove Mitochondrial, Ribosomal, Housekeeping, and lncRNA Genes using :py:func:`tardis_spac.utils.remove_mito_ribo_hk_lnc_genes()`

.. note::
   Housekeeping gene list is from `He2020Nature_mouseHK.txt <https://github.com/zenglab-pku/SPAC-seq/blob/main/He2020Nature_mouseHK.txt>`_.
   Refer to `Nature paper <https://github.com/brianpenghe/Matlab-genomics/blob/master/He_2020_ENCODE3_RNA/GeneLists/Bulk Cluster Ubiquitous.txt>`_ for more details.

.. py:function:: tardis_spac.utils.remove_mito_ribo_hk_lnc_genes(adata, housekeeping_list="He2020Nature_mouseHK.txt")

   Remove genes related to mitochondria, ribosomes, housekeeping, or annotated as lncRNAs from the dataset.

   :param adata: AnnData object. The input data matrix with genes in adata.var_names.
   :param housekeeping_list: Path to a text file containing housekeeping gene names. Default: "He2020Nature_mouseHK.txt"
   :return: AnnData object with unwanted genes removed.

Filter and Process Guide Reads using :py:func:`tardis_spac.utils.filter_guide_reads()`
This process is essential for removing noise guide reads and keeping major signals for more accurate analysis.

.. note::
   For Stereoseq data only.
   For Visium HD data, we recommend running **TARDIS** directly.

.. py:function:: tardis_spac.utils.filter_guide_reads(gem_path, guide_prefix=None, output_path=None, binarilize=False, assign_pattern='max', filter_threshold=None)

   Filter and process guide reads from a GEM file.

   :param gem_path: Path to the input GEM file.
   :param guide_prefix: Optional string. Filter guide names with this prefix.
   :param output_path: Path to save filtered results. If None, returns the filtered DataFrame.
   :param binarilize: If True, binarize all bin counts to 1 (recommended for high-resolution/high-library datasets).
   :param assign_pattern: How to handle bins with multiple guides: 'max' (keep guide with max count), 'drop' (remove such bins), or 'all' (keep all).
   :param filter_threshold: Integer. Minimum MIDCount threshold for bin filtering.
   :return: Filtered pandas DataFrame if output_path is None, otherwise None.

.. code-block:: python

   from tardis_spac.utils import filter_guide_reads
   filter_guide_reads(
      "guide.gem",
      guide_prefix="sg",
      output_path="filtered_guide.gem",
      binarilize=True, # Recommended for tumor / immune cell perturbation datasets
      assign_pattern="max", # Recommended for more clean guide assignment
      filter_threshold=10 # Minimum MIDCount (UMI) threshold for bin filtering
   )

Combine Guide Replicates using :py:func:`tardis_spac.utils.combine_guide_replicates()`

.. note::
   Guide replicates are guide probes targeting the same gene and are typically labeled with an underscore and a number.

   For example, a guide probe targeting the gene *A* with replicate number 1 would be labeled as *A_1*.

   Similarly, a guide probe targeting the gene *B* with replicate number 2 would be labeled as *B_2*.

   These guide probes are often used to measure the expression of the same gene in different conditions or tissues.

   In **TARDIS**, we collapse these guide probes into a single guide name by summing the counts of the replicate guides.

.. py:function:: tardis_spac.utils.combine_guide_replicates(gdata)

   Merge guide replicates in the input AnnData object by collapsing underscore-delimited replicates into a single guide name.

   :param gdata: AnnData object of guide count data with var_names containing replicates.
   :return: AnnData object with guides collapsed (replicates summed per guide).
