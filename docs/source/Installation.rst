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

To install TARDIS (R version), first install *devtools*

.. code-block::

   $ install.packages('devtools')
   $ devtools::install_github("zenglab-pku/TARDIS")

R will automatically install all the dependencies.

.. note::

   R version *4.2.2* is recommended.

.. [1] Varrone, M., Tavernari, D., Santamaria-Martínez, A. et al. CellCharter reveals spatial cell niches associated with tissue remodeling and cell plasticity. Nat Genet 56, 74–84 (2024).