=======================================================
riboraptor : Tool for analysing ribosome profiling data
=======================================================


.. image:: https://img.shields.io/pypi/v/riboraptor.svg
        :target: https://pypi.python.org/pypi/riboraptor

.. image:: https://travis-ci.com/saketkc/riboraptor.svg?token=GsuWFnsdqcXUSp8vzLip&branch=master
    :target: https://travis-ci.com/saketkc/riboraptor

.. image:: https://pyup.io/repos/github/saketkc/riboraptor/shield.svg
     :target: https://pyup.io/repos/github/saketkc/riboraptor/
     :alt: Updates


Python package to analyse ribosome profiling data


* Free software: BSD license
* Documentation: http://saketkc.github.io/riboraptor/

Installation
------------

Using conda:

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


====
TODO
====

riboraptor
----------


- [ ] Add Snakemake and config template
- [ ] Add tests with demo sra files
- [ ] Make conda requirement strict; Add more detailed installation docs
- [ ] Add annotation BED files for each genome's GTF and associated notebook
- [ ] Improve documentation:
    - [ ] QC - Fragment length, Periodicity, Summary of STAR mapping results
    - [ ] CDS/(UTR5+CDS+UTR3) distribution
    - [ ] Instructions for parsing gene_coverages.tsv file in both Python and R
- [ ] Add differential translation analysis with Riborex
- [ ] Add ribotapor for uORF or dORF identification

ribodb
------

- [ ] Genome browser track links organized by each project (or publication bibtex key)
- [ ] Visualization of fragment length distribution
- [ ] Visualization of metagene plots
- [ ] Visualization of genewise coverage
- [ ] Visualization of CDS/UTR enrichment across genes
- [ ] MultiQC report 

