=============================================================
riboraptor : a pipeline for analysing ribosome profiling data
=============================================================


.. image:: https://img.shields.io/pypi/v/riboraptor.svg?style=flat-square
        :target: https://pypi.python.org/pypi/riboraptor

.. image:: https://travis-ci.org/saketkc/riboraptor.svg?branch=master&style=flat-square
        :target: https://travis-ci.org/saketkc/riboraptor

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square
        :target: http://bioconda.github.io/recipes/riboraptor/README.html

.. .. image:: https://pyup.io/repos/github/saketkc/riboraptor/shield.svg
     :target: https://pyup.io/repos/github/saketkc/riboraptor/
     :alt: Updates

.. _Miniconda: https://conda.io/miniconda.html
.. _`aspera connect`: http://downloads.asperasoft.com/en/downloads/8?list
.. _`Line4 snakemake/jobscript.sh`: https://github.com/saketkc/riboraptor/blob/47c8a50753c2bcc96b57d43b525a47bb8fde2d04/snakemake/jobscript.sh#L4
.. _`Line6 snakemake/cluster.yaml`: https://github.com/saketkc/riboraptor/blob/47c8a50753c2bcc96b57d43b525a47bb8fde2d04/snakemake/cluster.yaml#L6
.. _`Line7 snakemake/cluster.yaml`: https://github.com/saketkc/riboraptor/blob/47c8a50753c2bcc96b57d43b525a47bb8fde2d04/snakemake/cluster.yaml#L7


Python package to analyse ribosome profiling data


* Free software: BSD license
* Documentation: https://riboraptor.readthedocs.io/en/latest/


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

   conda create --name riboraptor python=3.6 gcc matplotlib numpy pandas pybedtools \
   pyBigWig pyfaidx pysam scipy seaborn statsmodels six click click-help-colors htseq biopython \
   snakemake sra-tools star fastqc trim-galore ucsc-bedgraphtobigwig ucsc-bedsort r-rcurl \
   r-rsqlite r-devtools r-optparse bioconductor-biocinstaller bioconductor-annotationdbi \
   bioconductor-geometadb bioconductor-geoquery

Installing riboraptor
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   source activate riboraptor
   git clone git@github.com:saketkc/riboraptor.git
   cd riboraptor
   python setup.py install --single-version-externally-managed --record=record.txt


Using riboraptor
----------------

Usage mode 1: use riboraptor as a Snakemake based workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#. Create a copy of `snakemake/configs/template.py` say `snakemake/configs/myNewProject.py`
#. Edit the paths inside the config file to your RAW data, output, GTF and BED files
#. Export your miniconda path by editing `Line4 snakemake/jobscript.sh`_
#. Edit `Line6 snakemake/cluster.yaml`_ to your error log file
#. Edit `Line7 snakemake/cluster.yaml`_ to your output log file
#. Submit job

.. code-block:: bash

   cd snakemake
   bash submitall.sh myNewProject

Usage mode 2: use riboraptor as a standalone toolkit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See: https://riboraptor.readthedocs.io/en/latest/
  
Usage mode 3: use riboraptor in a Galaxy environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Check: http://nucleus.usc.edu:8080/


Usage mode 4: ribopod - database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In progress: http://nucleus.usc.edu:8050/

Downloading datasets from SRA
-----------------------------

#. Install `aspera connect`_ 

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

See: https://riboraptor.readthedocs.io/en/latest/cmd-manual.html


