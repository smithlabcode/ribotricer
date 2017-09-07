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


Features
--------

* TODO

Scores
------

* Ribosome Release Score [Guttman2013]: A score to determine if protein translation is complete.  Defined as the ratio between reads in coding region to reads in the 3'UTR region normalized by the corresponding ratio in mRNA data allowing for incorrect 3'UTR annotations (CDS to 3'UTR ratio in mRNA is not 1 and is reflective of different protein producing potential)

    $\frac{\frac{Counts_{CDS}}{Counts_{UTR}}_{Ribo}}{\frac{Counts_{CDS}}{Counts_{UTR}}_{RNA}}$


Credits
---------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. [Guttman2013]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3756563/
