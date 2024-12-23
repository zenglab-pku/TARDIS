Welcome to SPOT!
=======================================

This is the documentation for SPOT, a Python package for spatial perturbation analysis.

.. image:: https://img.shields.io/pypi/v/spot-seq.svg
   :target: https://pypi.python.org/pypi/spot-seq

.. image:: https://readthedocs.org/projects/spot-seq/badge/?version=latest
   :target: https://spot-seq.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. image:: https://img.shields.io/github/license/ZengLab/SPOT
   :target: https://github.com/ZengLab/SPOT/blob/main/LICENSE
   :alt: License

.. image:: https://img.shields.io/github/stars/ZengLab/SPOT?style=social
   :target: https://github.com/ZengLab/SPOT
   :alt: GitHub stars

Read the tutorial here:

https://spot-tutorial.readthedocs.io/en/latest/

Installing using python
------------

To use SPOT (python version), first install its dependencies using pip or conda.

You can simply install SPOT dependencies using pip through:

.. code-block:: 

   (.venv) $ pip install numpy pandas scipy scikit-learn anndata scanpy

**SPOT** can be ran on both Linux and Mac OS, since it is interpreted by python and R.

Installing using R
------------

To install SPOT (R version), first intall *dev-tools*

.. code-block::

   $ install.packages('devtools')
   $ devtools::install_github("pkuTrasond/SPOT")

R will automatically install all the dependencies.

For more detailed information, please refer to the tutorial.