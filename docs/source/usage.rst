Usage
=====

.. note::

   This here is an API version of all the functions applicable in **SPOT**,
   For more detailed and thorough reference please visit our :ref:`tutorial` site.

Preprocessing
----------------

Spatial perturbation can be highly arbitrary if we cannot perform valid
preprocessing and filtering of low quality guides and bins. Refer to *Paper name*
for our in house filtering method.

**SPOT** performs filtering with validation panels with the following methods.

First import **SPOT** from library.

.. code-block::

   import spot as sp
   sp.set_random_seed(42)

In this tutorial we use our in house spatial transcriptomics data.
This data incorporates a library of **68 guides**, and is sequenced on **BGI Stereo-seq** platform.

.. code-block::

   # perform quality check from BGI stereo-seq GEM output
   sp.preprocessing.filter_qc_bins('A04091E1.gem')

*output:*

.. image:: ../_images/filter_qc.png
   :align: center
   

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

