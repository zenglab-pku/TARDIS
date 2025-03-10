Installation
=====

.. _Installation:

Installing using python
------------

To use TARDIS (python version), first install its dependencies using pip or conda.

.. attention:: 
   To run **TARDIS** with cluster dependent modules, you must specify cluster fields.
   One recommended spatial clustering method is *Cellcharter* [1]_.

   Alternatively, we provide a Non-negative Matrix Factorization (NMF) model for simple analysis.

You can simply install TARDIS dependencies using pip through:

.. code-block:: 

   (.venv) $ pip install numpy pandas scipy scikit-learn anndata scanpy

**TARDIS** can be ran on both Linux and Mac OS, since it is interpreted by python and R.

.. note::

   Both python and R version requires *python* to run **Cluster dependent** analysis. since
   they require *Cellcharter* to perform clustering.

   If other clustering method is performed, follow the instructions in :ref:`tutorial` to perform
   downstream analysis.

Installing using R
------------

To install TARDIS (R version), first intall *dev-tools*

.. code-block::

   $ install.packages('devtools')
   $ devtools::install_github("zenglab-pku/TARDIS")

R will automatically install all the dependencies.

.. note::

   R version *4.2.2* is recommended.

.. [1] Varrone, M., Tavernari, D., Santamaria-Martínez, A. et al. CellCharter reveals spatial cell niches associated with tissue remodeling and cell plasticity. Nat Genet 56, 74–84 (2024).