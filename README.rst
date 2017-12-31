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
   $ conda install numpy scipy pandas matplotlib seaborn statsmodels pycwt mtspec htseq six pybigwig click pysam
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

