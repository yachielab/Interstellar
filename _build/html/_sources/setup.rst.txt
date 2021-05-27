Setup
-------

System requirements
~~~~~~~~~~~~~~~~~~

Barista was tested with following environment. Install packages using conda or pip.

- Python 3.7.3
- kneed 0.5.0
- matplotlib 3.0.3
- numpy 1.16.2
- pandas 0.24.2
- regex 2019.08.19
- symspellpy 6.5.2

Barista setup
~~~~~~~~~~~~~

Clone this repository:

.. code:: bash
    
    git clone /URL/to/barista.git

Then add the excutable authority to barista/script/barista and add the barista/script directory to $PATH:

.. code:: bash

    cd barista
    chmod u+x script/barista

Test run
~~~~~~~~

.. code:: bash

    barista -v
    >> barista version 1.0

