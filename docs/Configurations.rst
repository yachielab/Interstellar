Configurations
==============

Data processing configuration file
**********************************

The data processing configuration file consists of five setting sections, “general”, “value_extraction”, “value_translation”, 
“demultiplex” and “tag”. 

-------------------------------------

general (required)
^^^^^^^^^^^^^^^^^^

In this section, users configure general settings for command execution. 

Example
-------

.. code::

    [general]
    PROJECT_DIR=/path/to/project_directory
    SET_SHELL_ENV=/path/to/template.sh
    SAMPLESHEET=

Description
-----------

* **PROJECT_DIR (required):** Path to a project directory. All processes are performed under this directory.
* **SET_SHELL_ENV (required):** Path to a template shell script for your environment. Example file
* **SAMPLESHEET (optional):** A tsv file to specify the correspondence between the input file prefixes and their sample names. File prefix is defined as the unique part of the file name among FASTQ files for the same read. By providing this table, at the value translation step, Interstellar assigns distinct values to the same sequences in different samples. This function is useful when users want to translate multiple FASTQ files demultiplexed by sample indices and merge the emulated files into a single file. If there are three paired-end reads each of which corresponds to a distinct sample (plate10, plate3 and plate4) as follows:

  .. code::

      Read1
      ├── svf_plate10_SRR5664329_1.fastq.gz
      ├── svf_plate3_SRR5664331_1.fastq.gz
      └── svf_plate4_SRR5664330_1.fastq.gz

      Read2
      ├── svf_plate10_SRR5664329_2.fastq.gz
      ├── svf_plate3_SRR5664331_2.fastq.gz
      └── svf_plate4_SRR5664330_2.fastq.gz

  Then, the ``SAMPLESHEET`` should be specified like this:

  .. code::

      svf_plate10_SRR5664329	plate10
      svf_plate3_SRR5664331	plate3
      svf_plate4_SRR5664330	plate4

-----------------------------------------------

value_extraction (required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^
In this section, users configure how sequencing reads are interpreted. 
This section exports a corrected sequence segment table as a tsv format. 
Read processing with this section is required for the downstream processes, value_translation, demultiplex, and tag.
See the details.

Example
-------

.. code::

    [value_extraction]
    # Input file path
    READ1_DIR=/path/to/Interstellar/example-dataset/fastqs/R1
    READ2_DIR=/path/to/Interstellar/example-dataset/fastqs/R2

    # Read structure configuration
    READ1_STRUCTURE=^(?P<segment1>.{9,10})(?P<segment2>CAGAGC){s<=2}(?P<segment3>.{8})(?P<segment4>.{10}).*$
    READ2_STRUCTURE=(?P<segment5>.*)

    # Parental-local variable allocations
    segment1,segment3=CHILD_OF(segment1,segment4)

    # Quality filtering and correction for each sequence segment
    cb_value1=QUALITY_FILTER(source:segment1,min_nucleotide_Q-score:5,min_avg_Q-score:20) 
    >> WHITELIST_CORRECT(levenshtein_distance:1,path:/path/to/Intersellar/example-dataset/sci-RNA-seq3_HP_purified_whitelist.txt) 
    >> SEQ2VALUE()

    cb_value2=QUALITY_FILTER(source:segment4, min_nucleotide_Q-score:5, min_avg_Q-score:20) 
    >> KNEE_CORRECT(rank:auto, levenshtein_distance:1) 
    >> WHITELIST_CORRECT(levenshtein_distance:1, path:example-dataset//sci-RNA-seq3_RT_purified_whitelist.txt) 
    >> SEQ2VALUE()

    const_value4=PASS(source:segment2)

    umi_value3=QUALITY_FILTER(source:segment3, min_nucleotide_Q-score:5, min_avg_Q-score:20) >>  SEQ2VALUE()

    bio_value5=QUALITY_FILTER(source:segment5, min_nucleotide_Q-score:5, min_avg_Q-score:20)

Description
-----------

* **READ1_DIR, READ2_DIR, INDEX1_DIR, INDEX2_DIR (required at least one of them):** Path to directories containing input FASTQ files. FASTQ file names must follow this format:
    
  .. code::
    
    <fileprefix>_<read identifier>.fastq (fastq.gz is also allowed)
    <fileprefix> must be unique across the FASTQ files from the same read, and <read identifier> must be unique across 
    the FASTQ files from different reads with the same prefix. For example,

  * If the process is for multiple samples

    For multi-sample analysis, a ``SAMPLESHEET`` should be provided in the “general” section to identify the correspondence 
    between the files and the samples. The file prefixes in the ``SAMPLESHEET`` are used to collect the input files. 
    For the example ``SAMPLESHEET``, the input file directory should look like this:

    .. code::

        Input_dir
        ├── read1
        │   ├── svf_plate10_SRR5664329_1.fastq.gz
        │   ├── svf_plate3_SRR5664331_1.fastq.gz
        │   └── svf_plate4_SRR5664330_1.fastq.gz
        └── read2
            ├── svf_plate10_SRR5664329_2.fastq.gz
            ├── svf_plate3_SRR5664331_2.fastq.gz
            └── svf_plate4_SRR5664330_2.fastq.gz

  * If the process is for a single sample

    File prefixes are automatically collected from the input directories. The input directory should look like this:

    .. code::

        Input_dir/
        ├── read1
        │   ├── hgmm_1k_v3_S1_L001_R1_001.fastq.gz
        │   └── hgmm_1k_v3_S1_L002_R1_001.fastq.gz
        └── read2
            ├── hgmm_1k_v3_S1_L001_R2_001.fastq.gz
            └── hgmm_1k_v3_S1_L002_R2_001.fastq.gz


* **FLASH (optional):** Paired-end reads can be merged by specifying the paired reads here. The pair should be chosen from ``READ1-READ2``, ``READ1-INDEX1``, ``READ1-INDEX2``, ``READ2-INDEX1``, ``READ2-INDEX2``, or ``INDEX1-INDEX2``.
* **FLASH_MIN_OVERLAP (optional):** FLASH2 parameter corresponding to ``-m`` option to specify minimum overlap length. Default=20.
* **FLASH_MAX_OVERLAP (optional):** FLASH2 parameter corresponding to ``-M`` option to specify maximum overlap length. Default=30.
* **READ1_STRUCTURE:**  Read1 structure specified regular expression patterns. To set multiple patterns for a single read, add unique suffixes to the option names such as ``READ1_STRUCTURE_1``, ``READ1_STRUCTURE_2``, … Group names defined by ``?P<>`` are recognized as segment names. Interstellar adopts extended regular expression patterns implemented in `regex <https://pypi.org/project/regex/>`_ library, which allows mismatches in a constant sequence. Examples for the read structure definition are as follows:
  
  * ex) Quartz-seq library
 
    .. code::

        READ1_STRUCTURE=^(?P<seg1>.{14})(?P<seg2>.{8})$
        READ2_STRUCTURE=^(?P<seg3>.*)$

    In this example, read1 is 22 bp and it consists of 14 bps seg1 segment and 8 bp seg3 segment from 5’ to 3’. For read2, the entire sequence is recognized as seg3 segment. Here, seg1, seg2, and seg3 correspond to a cell barcode, a UMI, and a cDNA, respectively.

    .. image:: img/readme_regex.001.trim.jpeg
        :scale: 30 %

  * ex) Read1 of sci-RNA-seq3 library

    .. code::

        READ1_STRUCTURE=^(?P<seg1>.{9,10})(?P<seg2>CAGAGC){s<=2}(?P<seg3>.{8})(?P<seg4>.{10})$

    In this regular expression pattern, seg2 is first searched from a given sequence allowing a maximum of two mismatches. Then, seg1 is extracted if the upstream of the seg2 is 9 or 10 bp. Finally, 8 bp downstream of the seg2 and 10 bp downstream of the seg3 are extracted as seg3 and seg4, respectively. seg1, seg2, seg3, and seg4 respectively correspond to a cell barcode added by the ligation, a constant sequence, a UMI, and a cell barcode by the reverse transcription.

    .. image:: img/readme_regex.002.jpeg
        :scale: 30 %

* **READ2_STRUCTURE:** Read2 structure. See the disctiption of READ1_STRUCTURE.
* **INDEX1_STRUCTURE:** Index1 structure. See the disctiption of READ1_STRUCTURE.
* **INDEX2_STRUCTURE:** Index2 structure. See the disctiption of READ1_STRUCTURE.
* **READ_FLASH_STRUCTURE:** Read structure merged by FLASH2. See the disctiption of READ1_STRUCTURE.
* **Parental-local value allocations:** Using the segments defined in the regular expression patterns, users can define hierarchical relationships of the segments that we call ‘parental-local value allocations’ for optimal value conversion. In brief, under this allocation, local values have valid values in combination with its parental values on the same sequencing read. Therefore, the same local values can be assigned among the local values linked to distinct parental values, which can compress the information space. For example, unique molecular identifiers (UMIs) and cell barcodes in scRNA-seq correspond to local and parental values, respectively. UMIs encode values of distinct mRNA molecules in each cell, and two UMI values in different cells can theoretically be identical. In the configuration file, the allocation is specified by a function ``CHILD_OF()``.

  * For example, in case that ``seg1`` is a local value of ``seg2``:

      .. code::

          seg1 = CHILD_OF(seg2)
    
  * In case ``seg1`` is a local value of the combination of ``seg2`` and ``seg3`` (An example situation is UMI and cell barcodes in sci-RNA-seq – A combination of two cell barcode segments identifies a single cell):
    
      .. code::

          seg1 = CHILD_OF(seg2,seg3)

  * In case a combination of ``seg1`` and ``seg2`` is allocated under a combination of ``seg3`` and ``seg4``:

      .. code::

          seg1,seg2 = CHILD_OF(seg3,seg4)

* **Quality filtering, sequence error-correction, and transforming sequences to values:** Here, data processing strategies for each segment are described by functions of ``QUALITY_FILTER()``, ``KNEE_CORRECT()``, ``WHITELIST_CORRECT()``, ``BARTENDER()`` and ``SEQ2VALUE()``, and a pipe-like operator ``>>``.

  * ``QUALITY_FILTER()``
    
    This function performs sequence filtering based on their quality score. It accepts the following three arguments:

    * ``source`` (required): Target segment name.
    * ``min_nucleotide_Q-subdirectories``: Minimum allowable quality score across all bases in the segment. Default=5.
    * ``min_avg_Q-score``: Minimum allowable average quality score of the segment. Default=20.

  * ``KNEE_CORRECT()``
    
    This function performs sequence error correction based on the knee point in the rank-count plot of the segment. It accepts the following three arguments:

    * ``source``: Target segment name. If this argument is already specified at the upstream of >>, this is not required.
    * ``rank``: Method to determine the rank corresponding to the knee point in the rank-count plot. If ‘auto’ is specified, the knee point is automatically computed by kneed. A numerical value can also be provided to manually indicate the rank for the knee point. Default=auto.
    * ``levenshtein_distance``: Maximum Levenshtein distance to correct the minority sequences to the nearest sequences in the majority sequence population. Default=1.

  * ``WHITELIST_CORRECT()``
    
    This function performs sequence error correction based on the sequence whitelist. It accepts the following three arguments:

    * ``source``: Target segment name. If this argument is already specified at the upstream of >>, this is not required.
    * ``path`` (required): Path to the sequence whitelist.
    * ``levenshtein_dist``: Maximum Levenshtein distance to correct the sequences to the nearest sequences in the whitelist. Default=1.
 
  * ``BARTENDER()``
    
    This function performs sequence error correction using Bartender. It accepts the following six arguments:
    
    * ``source``: Target segment name. If this argument is already specified at the upstream of >>, this is not required.
    * ``-c``: The sequence frequency cutoff. See the details in the Bartender manual.
    * ``-d``: The maximum cluster distance that may be merged. See the details in the Bartender manual.
    * ``-z``: The cluster merging threshold. See the details in the Bartender manual.
    * ``-l``: The seed length. See the details in the Bartender manual.
    * ``-s``: The number of non-overlapping positions between two adjacent seeds. See the details in the Bartender manual.

  * ``SEQ2VALUE()``
    
    This function converts sequences into values for the target segment. This process is required if the segment is supposed to be converted in the value_translation step. It accepts the following argument:
    
    * ``source``: Target segment name. If this argument is already specified at the upstream of ``>>``, this is not required.
    * Pipe-like operator ``>>``: A pipe-like operator ``>>`` can connect multiple processes by the functions above and describe a process flow.

  Example of a process flow:

  .. code::

      segment1.filtered.corrected.value=QUALITY_FILTER(source:segment1,min_nucleotide_Q-score:5,min_avg_Q-score:20) 
      >> WHITELIST_CORRECT(levenshtein_distance:1,path:/path/to/whitelist.txt) 
      >>  SEQ2VALUE()

  In this example, the processed segment ‘segment1.filtered.corrected.value’ is derived from the segment ‘segment1’ followed by quality filtering, sequence error correction with the whitelist, and converting the sequences into values. This process flow can also be described redundantly as follows:

  .. code::

      segment1.filtered=QUALITY_FILTER(source:segment1,min_nucleotide_Q-score:5,min_avg_Q-score:20) 
      segment1.filtered.corrected=WHITELIST_CORRECT(source: segment1.filtered, levenshtein_distance:1,path:/path/to/whitelist.txt) 
      segment1.filtered.corrected.value=SEQ2VALUE(source:segment1.filtered.corrected)


---------------------------------------------






value_translation (optional)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this section, users configure the value translation strategies and the emulated read structure. 
Here, the source segment values are optimized considering the value combinations and the parental-local value allocations. 
Then, destination sequences are assigned to the optimized values, and emulated FASTQ files are exported. See the details.

.. _labelCap1:

Example
-------

.. code::

    [value_translation]

    #Value to destination sequence conversion
    dest_segment1=VALUE2SEQ(source:cb_value1,cb_value2, whitelist_path:/path/to/whitelist.txt)
    dest_segment2=VALUE2SEQ(source:umi_value3, randseq_length:12)
    dest_segment3=bio_value5
    exConst1=”TCTC”

    #Read structure configuration
    READ1_STRUCTURE=dest_segment1+dest_segment2+exConst1
    READ2_STRUCTURE=dest_segment3


Description
-----------

* **Defining destination segment sequences:** The destination segment sequences are defined by data processing functions, ``VALUE2SEQ()`` or ``SEQ2SEQ()``, or directory from a source segment. A single destination segment can be defined as:

  .. code::
    
    dest_segment = …

  If a combination of destination segments corresponds to source segments (The example situation is cell barcode conversion from 10X to SPLiT-seq), the destination segment can be defined by separating the destination segments by “,”:

  .. code::

    dest_segment1,dest_segment2 = …

  Then, at the read structure configuration step, you can use dest_segment1 and dest_segment2. This combinatorial destination segment is only applicable for ``VALUE2SEQ()`` and ``SEQ2SEQ()``.

  Examples for defining various destination sequences:

  ex1) one destination segment, three source segments, sequence assignments from a whitelist

    .. code::

      dest=VALUE2SEQ(source:seg1,seg2,seg3, path:/path/to/whitelist.txt)

  ex2) two destination segments, two source segments, random sequence assignments

  In this example, 4 bp and 6 bp random sequences are assigned to dest1 and dest2, respectively.

    .. code::

      dest1,dest2=VALUE2SEQ(source:seg1,seg2, randseq_length: 4,6)

  ex3) two destination segments, two source segments, conversion table-based assignments

    .. code::
    
      dest1,dest2=SEQ2SEQ(source:seg1,seg2, conversion_table: /path/to/table)

  * ``VALUE2SEQ()``

    This function assigns sequences from a given whitelist or random sequences to the destination segment. It accepts the following three arguments:

    * ``source`` (required): Target source segment name. To set the source as a combination of multiple source segments, specify those segments split by “,”.
    * ``whitelist_path``: Path to the whitelist to assign sequences to the destination segment. If multiple destination segments are defined, the same number and order of whitelist paths with the destination segments must be provided by splitting the paths by “,”.
    * ``randseq_length``: Random sequence length. If multiple destination segments are defined, the same number and order of the sequence lengths with the destination segments must be provided by splitting the lengths by “,”.

  * ``SEQ2SEQ()``

    This function assigns destination sequences from a user-defined sequence conversion table. An example situation is the spatial barcode conversion from Slide-seq to Visium where barcode correspondence must be preserved to preserve the coordinate information.

    * ``source`` (required): Target source segment name. To set the source as combination of multiple source segments, specify those segments split by “,”.
    * ``conversion_table``: Path to the tsv file that specifies the sequence correspondence between the source and the destination segment. Headers must correspond to the target source segments and the destination segments. In the following example table, the sequence combinations of source1 and source2 correspond to the combinations of dest1 and dest2.
    
    .. code::

      source1      source2     dest1         dest2
      ACTTGATTGT  TCCTACCAGT  ACTGAGTGC    CGAGATGT
      ACTTGATTGT  GCGTTGGAGC  TGATCGTCATG  GCGCATTC
      ACTTGATTGT  GATCTTACGC  GATCGGTACG   AACCCTTG
      ACTTGATTGT  CTGATGGTCA  GGCCCTTA     CGGACAAC
      ACTTGATTGT  CCGAGAATCC  ATCTGCATC    AGTTTAGA
      ACTTGATTGT  GCCGCAACGA  GATACAAACT   ATATACCT
      ACTTGATTGT  TGAGTCTGGC  TTATAGCC     GGTGGGAT
      ACTTGATTGT  TGCGGACCTA  AATTCAGTA    AAGGTAAT
      ACTTGATTGT  ACGGAGGCGG  GATCCGTCCA   ATGTTGGC

    

  * Constant sequence
    
    Constant sequence is specified by a nucleotide sequence sandwiched by single- or double-quotations (See :ref:`labelCap1` above).

  * Inheriting sequences from a source segment
    
    To inherit sequences from a source segment, users can specify the source segment name (See :ref:`labelCap1` above).

* **READ1_STRUCTURE:** Read1 structure. The structure should be descrived using the destination segments concatenated by “+”. 
  
* **READ2_STRUCTURE:** Read2 structure. See the disctiption of READ1_STRUCTURE.

* **INDEX1_STRUCTURE:** Index1 structure. See the disctiption of READ1_STRUCTURE.

* **INDEX2_STRUCTURE:** Index2 structure. See the disctiption of READ1_STRUCTURE.
  






------------------------------------

demultiplex (optional)
^^^^^^^^^^^^^^^^^^^^^^
In this section, users configure a strategy to demultiplex reads. 
Target sequence segments are grouped by user-specified segment sequences (key) and each group is separately exported 
in a FASTQ or tsv format. This method is performed after value extraction such that the demultiplexed reads or 
sequence segment table consist of quality-filtered and error-corrected segment sequences. See the details.

Example
-------

.. code::

    [demultiplex]
    KEY=seg1,seg2
    FORMAT=fastq
    TARGET=
    READ1_STRUCTURE=seg1+seg3
    READ2_STRUCTURE=seg4+seg5



Description
-----------

* **KEY (required):** Source segments to separate the file. If you need to split the files by combinations of multiple segments, specify them separated by “,”.

* **FORMAT (required):** Export data format selected from “FASTQ” or “tsv”.

* **TARGET:** Source segments to be included in the demultiplexed sequence table. This option is available if FORMAT=tsv.
  
* **READ1_STRUCTURE:** Read1 structure. The structure should be descrived using the destination segments concatenated by “+”. Not availbale if FORMAT=tsv.

* **READ2_STRUCTURE:** Read2 structure. See the disctiption of READ1_STRUCTURE.

* **INDEX1_STRUCTURE:** Index1 structure. See the disctiption of READ1_STRUCTURE.

* **INDEX2_STRUCTURE:** Index2 structure. See the disctiption of READ1_STRUCTURE.
  



---------------------------------

tag (optional)
^^^^^^^^^^^^^^^^
In this section, users configure a strategy to add user-defined segment sequences on header lines of the FASTQ file. 
This method is performed after value extraction such that the tagged reads consist of quality-filtered and error-corrected 
segment sequences. See the details.

Example
-------

.. code::

    [tag]
    READ1_TAG=seg1,seg2
    READ2_TAG=
    INDEX1_TAG=seg1,seg2,seg3
    INDEX2_TAG=
    READ1_STRUCTURE=seg1+seg3
    READ2_STRUCTURE=
    INDEX1_STRUCTURE=seg4+seg5
    INDEX2_STRUCTURE=



Description
-----------

* **READ1_TAG:** Source segments to be tagged on the output FASTQ headers of Read1. Multiple segments can be tagged by specifying them separated by “,”.

* **READ2_TAG:** Source segments to be tagged on the output FASTQ headers of Read2.

* **INDEX1_TAG:** Source segments to be tagged on the output FASTQ headers of Index1.

* **INDEX2_TAG:** Source segments to be tagged on the output FASTQ headers of Index2.

* **READ1_STRUCTURE:** Output Read1 structure. The structure should be descrived using the destination segments concatenated by “+”. Not availbale if FORMAT=tsv.

* **READ2_STRUCTURE:** Output Read2 structure. See the disctiption of READ1_STRUCTURE.

* **INDEX1_STRUCTURE:** Output Index1 structure. See the disctiption of READ1_STRUCTURE.

* **INDEX2_STRUCTURE:** Output Index2 structure. See the disctiption of READ1_STRUCTURE.

----------------------------------


Distributed computing configuration file
****************************************

This configuration file is used to perform distributed computing and is composed of one setting section, “qsub”. By providing this, Interstellar splits the input FASTQ files into multiple file splits and distributes the jobs described in the data processing configuration file in an SGE or UGE cluster environment.

Example
-------

.. code::

    [qsub]
    MEM_MAX=128
    MEM_RATIO=0.1
    QOPTION=-l s_vmem=<mem>G,mem_req=<mem>G
    NUM_READS=2000000


Description
-----------

* **MEM_MAX (required):** Maximum memory available in the environment. Numerical values should be provided.

* **MEM_RATIO:** A scaling factor for the memory requirement for the distributed jobs (0~1). This parameter is not obvious and depends on the read length, sequencing depth, and the number of segments. We recommend trying 0.05~0.1. The larger this parameter, the smaller numbers of distributed jobs are submitted, which will be inefficient. Default=0.05

* **QOPTION (required):** Options to be added to the qsub command. Resource specification options must be included since they depend on each cluster environment. The parts denoted by “<mem>” will be replaced by memory sizes.

* **NUM_READS:** A number of reads per split FASTQ file. A larger value decreases the number of the distributed jobs, but the computing time per job will be longer. Default=2000000.
