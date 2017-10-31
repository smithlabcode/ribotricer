===============================
Riboraptor
===============================


.. image:: https://img.shields.io/pypi/v/riboraptor.svg
        :target: https://pypi.python.org/pypi/riboraptor

.. image:: https://img.shields.io/travis/saketkc/riboraptor.svg
        :target: https://travis-ci.org/saketkc/riboraptor

.. image:: https://readthedocs.org/projects/riboraptor/badge/?version=latest
        :target: https://riboraptor.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/saketkc/riboraptor/shield.svg
     :target: https://pyup.io/repos/github/saketkc/riboraptor/
     :alt: Updates


Python package to analyse ribosome profiling data


* Free software: BSD license
* Documentation: https://riboraptor.readthedocs.io.

Installation
------------

.. code-block:: bash

   $ git clone git@github.com:saketkc/riboraptor.git
   $ conda install numpy scipy pandas matplotlib seaborn statsmodels pycwt
   $ cd riboraptor
   $ pip install -e .

Features
--------

* TODO

Scores
------

* **Ribosome Release Score** [Guttman2013]_ : 
.. image:: http://latex.codecogs.com/svg.latex?%5Cinline%20%5Cfrac%7B%5Cbig%28%5Cfrac%7BCounts_%7BCDS%7D%7D%7BCounts_%7BUTR%7D%7D%5Cbig%29_%7BRibo%7D%7D%7B%5Cbig%28%5Cfrac%7BCounts_%7BCDS%7D%7D%7BCounts_%7BUTR%7D%7D%5Cbig%29_%7BRNA%7D%7D
   :align: center
A score to determine if protein translation is complete.  Defined as the ratio between reads in coding region to reads in the 3'UTR region normalized by the corresponding ratio in mRNA data allowing for incorrect 3'UTR annotations (CDS to 3'UTR ratio in mRNA is not 1 and is reflective of different protein producing potential) 

* **ORFScore** [Bazzini2014]_ : Compares count of number of RPFs in each frame to a uniform distribution using Chi-Squared statistic to identify actively translated ORFs.
.. image:: https://latex.codecogs.com/svg.latex?%5Clog_2%5Cbig%281%20&plus;%20%5Csum_%7Bi%3D1%7D%5E3%20%5Cfrac%7B%28F_i-%5Cbar%7BF%7D%29%5E2%7D%7BF%7D%20%5Cbig%29%20%5Ctimes%20%5Cbegin%7Bcases%7D%20-1%20%26%20%28F_1%20%3C%20F2%29%20%5Ccup%20%28F_1%20%3C%20F_3%29%2C%5C%5C%201%20%26%20%5Ctext%7Botherwise%7D%20%5Cend%7Bcases%7D
where $F_i$ represents number of reads in frame $i$ and $\bar{F}$ represents $mean(F1,F2,F3)$

* **Floss Score** [Ingolia2014]_ : 
.. image:: https://latex.codecogs.com/svg.latex?0.5%20%5Ctimes%20%5Csum_%7Bl%3D26%7D%5E%7Bl%3D34%7D%20f%28l%29%20-%20f_%7Bref%7D%28l%29

where the $f_{ref}$ is contructed by counting the number of fragments of a particular read length over only annotated protein-coding genes. Cutiff is determined by identifying outliers using Tukey's method. 

    The FLOSS cutoff, calculated as a function of the total number of reads in the transcript histogram, was established by
    considering a rolling window of individual annotated genes and the computing the upper extreme outlier cutoff for each       window using Tukey’s method (Q3 + 3*IQR, where Q3 is the 3rd quartile and IQR is the interquartile range).


Credits
---------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. [Guttman2013] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3756563/
.. [Bazzini2014] https://www.ncbi.nlm.nih.gov/pubmed/24705786
.. [Ingolia2014] https://www.ncbi.nlm.nih.gov/pubmed/25159147
