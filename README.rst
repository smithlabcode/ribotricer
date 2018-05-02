=============================================================
riboraptor : a pipeline for analysing ribosome profiling data
=============================================================


.. image:: https://img.shields.io/pypi/v/riboraptor.svg?style=flat-square
        :target: https://pypi.python.org/pypi/riboraptor

.. image:: https://travis-ci.org/saketkc/riboraptor.svg?branch=master&style=flat-square
        :target: https://travis-ci.org/saketkc/riboraptor

.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square
        :target: http://bioconda.github.io/recipes/riboraptor/README.html

.. image:: https://readthedocs.org/projects/riboraptor/badge/?version=latest
        :target: http://riboraptor.readthedocs.io/en/latest/?badge=latest&style=flat-square


.. _Miniconda: https://conda.io/miniconda.html
.. _`aspera connect`: http://downloads.asperasoft.com/connect2/
.. _`Line4 snakemake/jobscript.sh`: https://github.com/saketkc/riboraptor/blob/47c8a50753c2bcc96b57d43b525a47bb8fde2d04/snakemake/jobscript.sh#L4
.. _`Line6 snakemake/cluster.yaml`: https://github.com/saketkc/riboraptor/blob/47c8a50753c2bcc96b57d43b525a47bb8fde2d04/snakemake/cluster.yaml#L6
.. _`Line7 snakemake/cluster.yaml`: https://github.com/saketkc/riboraptor/blob/47c8a50753c2bcc96b57d43b525a47bb8fde2d04/snakemake/cluster.yaml#L7
.. _`SRAdb`: https://www.bioconductor.org/packages/3.7/bioc/html/SRAdb.html
.. _`GSE37744`: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37744



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
We will create a searate environment inside conda for running `riboraptor`. The environment name is also `riboraptor`.
If you already have a conda environment named `riboraptor`, you can delete it by running:

.. code-block:: bash

   source deactivate riboraptor && conda env remove -n riboraptor

We will now install all the dependencies:

.. code-block:: bash

   conda create --name riboraptor python=3.6 gcc matplotlib numpy pandas pybedtools \
   pyBigWig pyfaidx pysam scipy seaborn statsmodels six click click-help-colors htseq biopython \
   trackhub snakemake sra-tools star fastqc trim-galore ucsc-bedgraphtobigwig ucsc-bedsort r-rcurl \
   r-rsqlite r-devtools r-optparse bioconductor-biocinstaller bioconductor-annotationdbi \
   bioconductor-geometadb bioconductor-geoquery && source activate riboraptor
  
We also have the following two dependencies for processing and downloading SRA datasets:
   
#. `aspera connect`_ : For allowing '.fasp' downloads from SRA

#. `SRAdb`_ : For fetching all experiments of a SRA project with the associated metadata

Since there is currently a bug in bioconductor-sradb, we will install it from github. 

.. code-block:: bash
   
   git clone https://github.com/seandavi/SRAdb
   cd SRAdb
   
Run `R`, and install SRAdb within `R` use `devtools`. Please make sure your `riboraptor` environment is already activated. (`source activate riboraptor`):

.. code-block:: r

   library(devtools)
   devtools::install(".")

And finally, we need two metadata files for processing SRA records:

.. code-block:: bash
    
   mkdir riboraptor-data && cd riboraptor-data
   wget -c http://starbuck1.s3.amazonaws.com/sradb/GEOmetadb.sqlite.gz && gunzip GEOmetadb.sqlite.gz
   wget -c https://starbuck1.s3.amazonaws.com/sradb/SRAmetadb.sqlite.gz && gunzip SRAmetadb.sqlite.gz
  

Installing riboraptor
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   source activate riboraptor
   git clone git@github.com:saketkc/riboraptor.git
   cd riboraptor
   python setup.py install --single-version-externally-managed --record=record.txt

We will assume you have the following directory structure for the rest of our analysis:

::

    | some_root_directory
    | ├── riboraptor
    | │   ├── snakemake
    | │   └── setup.py
    | ├── riboraptor-data
    | │   ├── GEOmetadb.sqlite
    | │   └── SRAmetadb.sqlite


Using riboraptor
----------------

Usage mode 1: use riboraptor as a Snakemake based workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



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


Features
--------

See: https://riboraptor.readthedocs.io/en/latest/cmd-manual.html


Example Workflow
----------------

We will be working with the first published Ribo-seq dataset `GSE37744`_ from Ingolia et al. (2012) which has three samples from JEK293 cell line profiled under different concentrations of magnesium (which affects the degree of digestion).

At this point, we assume you have already completed all the steps under `Installing dependencies`_. 


#. Step 1: Downloading datasets

We will download all SRA files corresponding to GSE13750.

.. code-block:: bash
   
   cd riboraptor
   download_sra_data --sradb=../riboraptor-data/SRAmetadb.sqlite \
   --geodb=../riboraptor-data/GEOmetadb.sqlite GSE37744

GEO IDs are automatatiicaly converted to corresponding SRP IDs. GSE37744 corresponds to SRP012648.
We will now use Snakemake to run all the downstream steps.

#. Step 2: Copy template

.. code-block:: bash
   
   cd snakemake
   cp configs/template.py.sample configs/SRP012648.py
   
   Edit the paths inside `SRP012648.py` to point to your RAW data, GTF and BED files.
   An example of a config would be:
   
.. code-block:: python
   
   ## Rawdata directory (as created by running download_sra_data)
   RAWDATA_DIR = '/staging/as/skchoudh/SRA_datasets/SRP012648'
   
   ## Output directory (will be created if does not exist)
   OUT_DIR = '/staging/as/skchoudh/riboraptor-analysis/SRP012648'   
   
   ## Genome fasta location
   GENOME_FASTA = '/home/cmb-06/as/skchoudh/genomes/hg38/fasta/hg38.fa'
   
   ## Chromosome sizes location
   CHROM_SIZES = '/home/cmb-06/as/skchoudh/genomes/hg38/fasta/hg38.chrom.sizes'
   
   ## Path to STAR index (will be generated if does not exist)
   STAR_INDEX = '/home/cmb-06/as/skchoudh/genomes/hg38/star_annotated'

   ## GTF path
   GTF = '/home/cmb-06/as/skchoudh/genomes/hg38/annotation/gencode.v25.annotation.without_rRNA_tRNA.gtf'

   ## Path to bed file containing CDS coordinates coordinates
   CDS_BED = '/home/cmb-panasas2/skchoudh/riboraptor/riboraptor/annotation/hg38/cds.bed'
   
   ## Path to bed file containing 5'UTR coordinates coordinates
   UTR5_BED = '/home/cmb-panasas2/skchoudh/riboraptor/riboraptor/annotation/hg38/utr5.bed'

   ## Path to bed file containing 3'UTR coordinates coordinates
   UTR3_BED = '/home/cmb-panasas2/skchoudh/riboraptor/riboraptor/annotation/hg38/utr3.bed'

#. Step 3 : Change your miniconda path in `Line4 snakemake/jobscript.sh`_

   An example path would be:
.. code-block:: bash  
    
   export PATH="/home/cmb-panasas2/wenzhenl/miniconda3/bin:$PATH"


#. Step 4: Edit `Line6 snakemake/cluster.yaml`_ and `Line7 snakemake/cluster.yaml`_ to point to your log directory error log file
  
   An example path would be:
   
.. code-block:: yaml

   logout: '/home/cmb-06/as/wenzhenl/logs/{rule}.{wildcards.sample}.out'
   logerror: '/home/cmb-06/as/wenzhenl/logs/{rule}.{wildcards.sample}.err'

You would want to just edit the directory path leading to `/home/cmb-06/as/wenzhenl/logs/` and leave the rest as it is.

#. Step 5: Submit job

.. code-block:: bash

   bash submitall.sh SRP012648

The `submitall.sh` looks for a file named `SRP012648.py` in the configs directory, so make sure `SRP012648.py` exists inside
`configs/` directory.


