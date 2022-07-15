.. _header-n0:

Installation
============

Requirements
************

Install the following libraries:

- `matplotlib <https://anaconda.org/conda-forge/matplotlib>`_
- `numpy <https://anaconda.org/conda-forge/numpy>`_
- `pandas <https://anaconda.org/conda-forge/pandas>`_
- `regex <https://pypi.org/project/regex/>`_
- `kneed <https://pypi.org/project/kneed/>`_
- `symspellpy <https://pypi.org/project/symspellpy/>`_

.. note::

  The following versions of libraries were used for testing Interstellar.

  - Python 3.7.3
  - matplotlib 3.4.2
  - numpy 1.20.2
  - pandas 1.2.4
  - regex 2021.4.4
  - kneed 0.7.0
  - symspellpy 6.7.0

Also, the following command-line tools are required. Add the paths to $PATH.

- `seqkit <https://bioinf.shenwei.me/seqkit/>`_
- `FLASH2 <https://github.com/dstreett/FLASH2>`_
- `bartender <https://github.com/LaoZZZZZ/bartender-1.1>`_

Interstellar setup
*******************


Clone the repository::

  git clone /URL/to/Interstellar.git


Then add the executable authority to Interstellar/py and add the Interstellar/py directory to $PATH::

  cd Interstellar
  chmod u+x py/*

Then, confirm if all the libraries are successfully imported::

  Interstellar -h


Distributed computing
**********************

To enable distributed computing, qsub and qacct commands must be available in the environment to submit and monitor the Interstellar jobs.
This mode was tested with UGE 8.6.12 and SGE 6.2u6.

