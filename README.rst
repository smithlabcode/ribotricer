=======================================================
riboraptor : a pipeline for analysing ribosome profiling data
=======================================================


.. image:: https://img.shields.io/pypi/v/riboraptor.svg
        :target: https://pypi.python.org/pypi/riboraptor

.. image:: https://travis-ci.org/saketkc/riboraptor.svg?branch=master
        :target: https://travis-ci.org/saketkc/riboraptor

.. image:: https://pyup.io/repos/github/saketkc/riboraptor/shield.svg
     :target: https://pyup.io/repos/github/saketkc/riboraptor/
     :alt: Updates

.. _Miniconda: https://conda.io/miniconda.html

Python package to analyse ribosome profiling data


* Free software: BSD license
* Documentation: http://saketkc.github.io/riboraptor/

Installation
============

Setting up conda
****************

#. Install `conda`, the best way to install it is with the Miniconda_ package.The Python 3 version is recommended.

#. Set up channels, **It is important to add them in this order**.

.. code-block:: bash

   conda config --add channels r
   conda config --add channels defaults
   conda config --add channels conda-forge
   conda config --add channels bioconda

Installing depende
**************************

#. Create Python 3 conda env with name (e.g. riboraptor)
.. code-block:: bash

   conda create --name riboraptor python=3.6

.. code-block:: bash

   $ git clone git@github.com:saketkc/riboraptor.git
   $ conda create --name raptor gffutils matplotlib mne mtspec numpy pandas pybedtools \
   pyBigWig pyfaidx pysam scipy seaborn pycwt six click click-help-colors htseq biopython
   $ cd riboraptor
   $ python setup.py install --single-version-externally-managed --record=/tmp/record.txt


Using pip only:

.. code-block:: bash

   $ git clone git@github.com:saketkc/riboraptor.git
   $ cd riboraptor
   $ pip install -r requirements.txt
   $ python setup.py install


Features
--------

See: http://saketkc.github.io/riboraptor/cmd-manual.html


