Installation
=====

.. _installation:

Installing using python
------------

To use SPOT (python version), first install it using pip:

.. code-block:: 

   (.venv) $ pip install 

**SPOT** can be ran on both Linux and Mac OS, since it is interpreted by python and R.

.. note::

   Both python and R version requires *python* to run **Cluster dependent** analysis. since
   they require *Cellcharter* to perform clustering.

   If other clustering method is performed, follow the instructions in :ref:`tutorial` to perform
   downstream analysis.

Installing using R
------------

To install SPOT (R version), first intall *dev-tools*

.. code-block::

   $ install.packages('devtools')
   $ devtools::install_github("pkuTrasond/SPOT")

.. note::

   R version *4.2.2* is reconmended.