Installation
=====

.. _Installation:

Installing using python
-----------------------

To use TARDIS (python version), first install its dependencies using pip or conda.

.. attention:: 
   To run **TARDIS** with cluster dependent modules, you must specify cluster fields.
   One recommended spatial clustering method is *Cellcharter* [1]_.

   Alternatively, we provide a Non-negative Matrix Factorization (NMF) model for simple analysis.

You can simply install TARDIS dependencies using conda through:

.. code-block:: 

   (.venv) $ conda create -n tardis_env python>=3.8
   (.venv) $ conda activate tardis_env
   (.venv) $ conda install numpy pandas scipy scikit-learn anndata scanpy
   (.venv) $ pip install tardis_spac

.. note::

   Both python and R version requires *python* to run **Cluster dependent** analysis. since
   they require *Cellcharter* to perform clustering.
   
   If other clustering method is performed, follow the instructions in :ref:`tutorial` to perform
   downstream analysis.

Installing using R
------------

To install TARDIS (R version), first download the package from `GitHub <https://github.com/stereoseq/TARDIS>`_

.. note::
   Compatible with **R >= 4.2**. Uses **Seurat** for count matrices and guide-level metadata.

   The Python **CellCharter** module is not included in this port.

Install with devtools
~~~~~~~~~~~~~~~~~~~~~

Install core dependencies first (recommended), then install the package:

.. code-block:: r

   install.packages(c("devtools", "Matrix", "ggplot2"))
   # Seurat: install from CRAN or https://satijalab.org/seurat/articles/install.html
   install.packages("Seurat")

   devtools::install("/path/to/stereoseq/tardis_r")
   # If dependency download fails, use:
   # devtools::install("/path/to/stereoseq/tardis_r", dependencies = FALSE)

Alternative without devtools:

.. code-block:: r

   install.packages("/path/to/stereoseq/tardis_r", repos = NULL, type = "source")

For **h5ad** files (e.g. ``filtered_guide_bc_matrix.h5``), see :ref:`_Tutorial` for more details. Also set Python with ``anndata``:

.. code-block:: r

   install.packages("reticulate")
   Sys.setenv(RETICULATE_PYTHON = "/path/to/python/with/anndata")
   reticulate::use_python(Sys.getenv("RETICULATE_PYTHON"), required = TRUE)

Optional: ``NMF``, ``dbscan``, ``ks``, ``transport``, ``hdf5r``.

Quick start
~~~~~~~~~~~

.. code-block:: r

   library(tardisSpac)

   guide <- read_guide_h5("filtered_guide_bc_matrix.h5")

   guide <- filter_guide_reads_h5(
     guide,
     guide_prefix = "sg",
     filter_threshold = 10,
     assign_pattern = "max"
   )

   guide <- kl_divergence(guide, reference_guide = "sgTgfbr2_1", copy = TRUE)
   plot_top_kde(guide, result_field = "kl_div")

Python ↔ R API
~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1

   * - Python (``tardis_spac``)
     - R (``tardisSpac``)
   * - ``AnnData``
     - Seurat object
   * - ``.var[metric]``
     - ``object[["RNA"]]@meta.features[[metric]]``
   * - ``.obs``
     - ``object@meta.data``
   * - ``obsm['spatial']``
     - ``object@misc$spatial``
   * - ``read_guide_h5`` / ``sc.read_h5ad``
     - ``read_guide_h5()``
   * - ``stats.kl_divergence``
     - ``kl_divergence()``
   * - ``stats.wasserstein_distance``
     - ``wasserstein_distance()``
   * - ``stats.aitchison_distance``
     - ``aitchison_distance()``
   * - ``stats.permanova``
     - ``permanova()``
   * - ``utils.filter_guide_reads``
     - ``filter_guide_reads()``
   * - ``utils.filter_guide_reads_h5``
     - ``filter_guide_reads_h5()``
   * - ``utils.load_bin``
     - ``load_bin()``
   * - ``external.cluster_cellcharter``
     - *(not ported)*

.. [1] Varrone, M., Tavernari, D., Santamaria-Martínez, A. et al. CellCharter reveals spatial cell niches associated with tissue remodeling and cell plasticity. Nat Genet 56, 74–84 (2024).