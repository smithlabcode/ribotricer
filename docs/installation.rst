.. highlight:: shell

============
Installation
============

Setting up conda
----------------

#. Install `conda`, the best way to install it is with the Miniconda_ package.The Python 3 version is recommended.

#. Set up channels, **It is important to add them in this order**.

.. code-block:: bash

   conda config --add channels r
   conda config --add channels defaults
   conda config --add channels conda-forge
   conda config --add channels bioconda

Installing dependencies
-----------------------

.. code-block:: bash

   conda create --name riboraptor python=3.6 matplotlib numpy pandas pybedtools \
   pyBigWig pyfaidx pysam scipy seaborn statsmodels six click click-help-colors htseq biopython \
   snakemake sra-tools star fastqc trim-galore ucsc-bedgraphtobigwig ucsc-bedsort

Installing riboraptor
---------------------

.. code-block:: bash

   source activate riboraptor
   git clone git@github.com:saketkc/riboraptor.git
   cd riboraptor
   python setup.py install --single-version-externally-managed --record=/tmp/record.txt

Stable release
--------------

To install riboraptor, run this command in your terminal:

.. code-block:: console

    $ pip install riboraptor

This is the preferred method to install riboraptor, as it will always install the most recent stable release. 

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


From sources
------------

The sources for riboraptor can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/saketkc/riboraptor

Or download the `tarball`_:

.. code-block:: console

    $ curl  -OL https://github.com/saketkc/riboraptor/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install


.. _Github repo: https://github.com/saketkc/riboraptor
.. _tarball: https://github.com/saketkc/riboraptor/tarball/master
