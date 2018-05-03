Welcome to riboraptor's documentation!
======================================

``riboraptor`` library is a Python library for analysis of Ribo-seq data.
It assumes the reads have already been aligned to a reference and are available as BAM/SAM for input.

Current capabilities of ``riboraptor`` include:

:Visualization:
    - Plot read length distribution
    - Plot metagene coverage

:Utilities:
    - Export all gene coverages
    - Export metagene coverage
    - Export read length distribution
    - Calculate periodicity

.. code-block:: console

    $ Usage: riboraptor [OPTIONS] COMMAND [ARGS]...

      riboraptor: Tool for ribosome profiling analysis

     Options:
        --version  Show the version and exit.
        --help     Show this message and exit.

     Commands:
        bam-to-bedgraph              Convert bam to bedgraph
        bedgraph-to-bigwig           Convert bedgraph to bigwig
        export-bed-fasta             Export gene level fasta from specified bed
        export-gene-coverages        Export gene level coverage for all genes for given region
        export-metagene-coverage     Export metagene coverage for given region
        export-read-length           Calculate read length distribution
        periodicity                  Calculate periodicity
        plot-metagene                Plot metagene profile
        plot-read-length             Plot read length distribution
        uniq-bam                     Create a new bam with unique mapping reads only
        uniq-mapping-count           Count number of unique mapping reads


Contents:

.. toctree::
   :maxdepth: 2

   installation
   example_workflow
   cmd-manual
   api-usage
   scores
   modules
   history
   contributing
   authors



===============================================================================

Site map
========
 - :ref:`genindex`
 - :ref:`modindex`
 - :ref:`search`

