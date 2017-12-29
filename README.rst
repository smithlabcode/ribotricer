===============================
Ribocop
===============================


.. image:: https://img.shields.io/pypi/v/ribocop.svg
        :target: https://pypi.python.org/pypi/ribocop

.. image:: https://travis-ci.com/saketkc/ribocop.svg?token=GsuWFnsdqcXUSp8vzLip&branch=master
    :target: https://travis-ci.com/saketkc/ribocop

.. image:: https://pyup.io/repos/github/saketkc/ribocop/shield.svg
     :target: https://pyup.io/repos/github/saketkc/ribocop/
     :alt: Updates


Python package to analyse ribosome profiling data


* Free software: BSD license
* Documentation: http://saketkc.github.io/ribocop/

Installation
------------

Using conda:
.. code-block:: bash

   $ git clone git@github.com:saketkc/ribocop.git
   $ conda install numpy scipy pandas matplotlib seaborn statsmodels pycwt mtspec htseq six pybigwig click pysam
   $ cd ribocop
   $ python setup.py install --single-version-externally-managed --record=/tmp/record.txt


Using pip only:
.. code-block:: bash

   $ git clone git@github.com:saketkc/ribocop.git
   $ cd ribocop
   $ pip install -r requirements.txt
   $ python setup.py install


Features
--------

See: http://saketkc.github.io/ribocop/cmd-manual.html

