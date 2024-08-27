Usage
=====

.. _installation:

Installation
------------

To use SPOT (python version), first install it using pip:

.. code-block:: 

   (.venv) $ pip install 

To install SPOT (R version), first intall *dev-tools*

.. code-block::

   $ install.packages('dev-tools')

**SPOT** can be ran on both Linux and Mac OS, since it is interpreted by python and R.

.. note::

   Both python and R version requires *python* to run **Cluster dependent** analysis. since
   they require *Cellcharter* to perform clustering.

   If other clustering method is performed, follow the instructions in :ref:`tutorial` to perform
   downstream analysis.

Creating recipes
----------------

.. image:: ../_images/filter_qc.png

To retrieve a list of random ingredients,
you can use the ``lumache.get_random_ingredients()`` function:

.. autofunction:: lumache.get_random_ingredients

The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
will raise an exception.

.. autoexception:: lumache.InvalidKindError

For example:

>>> import lumache
>>> lumache.get_random_ingredients()
['shells', 'gorgonzola', 'parsley']

