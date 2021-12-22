Execution
=========

General usages
**************

Interstellar can take two different configuration files, one specifies data processing settings provided by ``-conf`` (required) 
and the other does distributed computing settings provided by ``-qconf`` (optional). 
When a distributed computing configuration file is provided, Interstellar splits input FASTQ files into small file splits and 
performs data processing for each split by generating qsub jobs automatically.

To run Interstellar with a single computing node::

    Interstellar -conf process.conf

To run Interstellar with distributed computing:: 

    Interstellar -conf process.conf -qconf distributed_computing.conf

The grammer of the configure file contents are explained here.


Outputs
***************

All Interstellar processes are conducted in a user-defined project directory. 
Output files produced by each method are exported in ``out`` subdirectory in the corresponding method directory. 
This is an example project directory structure for single sample analysis:

.. code:: 

    Project_directory
    ├── _flg
    └── sample1
        ├── config
        ├── filesplit
        ├── qlog
        ├── sh
        ├── demultiplex
        │   ├── _work
        │   └── out
        ├── tag
        │   ├── _work
        │   └── out
        ├── value_extraction
        │   ├── _work
        │   └── out
        └── value_translation
            ├── _work
            └── out


Interstellar can perform four methods, value_extraction, value_translation, demultiplex, and tag, 
and each of the methods respectively exports the following files in the ``out`` subdirectories:

- value_extraction: A quality-filtered and error-corrected sequence segment table. 
- value_translation: Emulated FASTQ files.
- demultiplex: FASTQ or tsv files split by the user-defined segment sequences.
- tag: FASTQ files with user-defined segment sequences tagged on header lines
