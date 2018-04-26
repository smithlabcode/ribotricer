=============================================================
riboraptor : a pipeline for analysing ribosome profiling data
=============================================================


.. image:: https://img.shields.io/pypi/v/riboraptor.svg
        :target: https://pypi.python.org/pypi/riboraptor

.. image:: https://travis-ci.org/saketkc/riboraptor.svg?branch=master
        :target: https://travis-ci.org/saketkc/riboraptor

.. .. image:: https://pyup.io/repos/github/saketkc/riboraptor/shield.svg
     :target: https://pyup.io/repos/github/saketkc/riboraptor/
     :alt: Updates

.. _Miniconda: https://conda.io/miniconda.html
.. _`aspera connect`: http://downloads.asperasoft.com/en/downloads/8?list


Python package to analyse ribosome profiling data


* Free software: BSD license
* Documentation: http://saketkc.github.io/riboraptor/


Installation
------------

Setting up conda
~~~~~~~~~~~~~~~~

#. Install `conda`, the best way to install it is with the Miniconda_ package.The Python 3 version is recommended.

#. Set up channels, **It is important to add them in this order**.

.. code-block:: bash

   conda config --add channels r
   conda config --add channels defaults
   conda config --add channels conda-forge
   conda config --add channels bioconda

Installing dependencies
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   conda create --name riboraptor python=3.6 matplotlib numpy pandas pybedtools \
   pyBigWig pyfaidx pysam scipy seaborn statsmodels six click click-help-colors htseq biopython \
   snakemake sra-tools star fastqc trim-galore ucsc-bedgraphtobigwig ucsc-bedsort

Installing riboraptor
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   source activate riboraptor
   git clone git@github.com:saketkc/riboraptor.git
   cd riboraptor
   python setup.py install --single-version-externally-managed --record=/tmp/record.txt

Downloading datasets from SRA
-----------------------------

#. Install `aspera connect`_ 
#. Install additional dependencies

.. code-block:: bash

   source activate riboraptor
   conda install gcc
   conda install -c r r=3.4.1
   conda install -c bioconda bioconductor-annotationdbi bioconductor-geometadb
   conda install -c r r-devtools
 
Since there is currently bug with bioconductor-sradb, we will install it from github

.. code-block:: bash

   git clone https://github.com/seandavi/SRAdb
   cd SRAdb
   
Run `R`, and install SRAdb within `R` use `devtools`

.. code-block:: r

   library(devtools)
   devtools::install(".")


Features
--------

See: http://saketkc.github.io/riboraptor/cmd-manual.html


