.. _header-n0:

Import
======

``import`` command is used to extract sequnce fragments by specified
format. The read structure is parsed as regular expression as default,
but it can also be described as 'simple' format.

.. _header-n6:

Usage
-----

.. code:: shell

   barista \
   -cmd import \
   -conf config.conf \
   -R1 input_R1.fastq.gz \
   -R2 input_R2.fastq.gz \
   -fq_gz \
   -d output_directory \
   -o output_filename

+-----------------------+-----------------------+-----------------------+
|  Option               | Required              | Description           |
+=======================+=======================+=======================+
| -cmd                  | ○                     | barista command       |
|                       |                       | (import).             |
+-----------------------+-----------------------+-----------------------+
| -conf                 | ○                     | Setting file. See     |
|                       |                       | below how to write    |
|                       |                       | the config file.      |
+-----------------------+-----------------------+-----------------------+
| -R1                   | ×                     | Input fastq file,     |
|                       |                       | read1.                |
+-----------------------+-----------------------+-----------------------+
| -R2                   | ×                     | Input fastq file,     |
|                       |                       | read2.                |
+-----------------------+-----------------------+-----------------------+
| -I1                   | ×                     | Input fastq file,     |
|                       |                       | index1.               |
+-----------------------+-----------------------+-----------------------+
| -I2                   | ×                     | Input fastq file,     |
|                       |                       | index2.               |
+-----------------------+-----------------------+-----------------------+
| -fq_gz                | ×                     | The Input fastq files |
|                       |                       | are copressed as gzip |
|                       |                       | format. As default,   |
|                       |                       | fastq files are       |
|                       |                       | regarded as           |
|                       |                       | un-compressed.        |
+-----------------------+-----------------------+-----------------------+
| -simple               | ×                     | Read structure is     |
|                       |                       | described as 'simple' |
|                       |                       | format, not regular   |
|                       |                       | expression.           |
+-----------------------+-----------------------+-----------------------+
| -d                    | ×                     | Output directory.     |
+-----------------------+-----------------------+-----------------------+
| -o                    | ×                     | Output file name.     |
+-----------------------+-----------------------+-----------------------+

.. _header-n279:

Setting file
------------

.. code:: 

   [import]
   src_raw_components=segment1, segment2, segment3, segment4, segment5
   barcodes=segment1,segment4,segment3
   Read1_src=^(?P<segment1>.{9,10})(?P<segment2>CAGAGC){s<=2}(?P<segment3>.{8})(?P<segment4>.{10}).*$
   Read2_src=(?P<segment5>.*)
   Index1_src=
   Index2_src=

.. _header-n178:

Description
~~~~~~~~~~~

-  src\_raw\_components **(Required)**

   The name of each component used in this ``import`` section, listed in
   comma-separated format. Here, segment1~5 were set as the components
   in this sequencing platform (specicically, this is a sciRNA-seq3
   structure and seqment1~5 respectively corresponds to cell barcode1,
   constant sequence, UMI, cell barcode2 and biological sequence).

-  barcodes **(Required)**

   The names of components listed in comma-separated format whose
   sequences are to be corrected and converted into other barcode
   sequneces of the destination sequence structure in following steps.

-  Read1\_src, Read2\_src, Index1\_src, Index2\_src

   The sequence structure of each read written in regular expression.
   Use component names specified in ``src_raw_components`` for the
   corresponding groups. Here, first 9 or 10 bases (in this case, cell
   barcode1) are extracted as segment1, subsequent CAGAGC (constant) is
   extracted as segment2 allowing 2 mismatches and next 8 and 10 bases
   (UMI and cell barcode2) are respectivery extracted as segment 4 and 5
   in read1. In read2, whole sequences are assigned to
   segment5(biological sequence). See the details about fuzzy matching
   `here <https://pypi.org/project/regex/>`__.

.. _header-n217:

Note
----

.. _header-n223:

Multiple read structure patterns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

User can specify multiple patterns for the read structures.

.. code:: shell

   [import]
   src_raw_components=segment1, segment2, segment3, segment4, segment5
   barcodes=segment1,segment4,segment3
   Read1_src=^(?P<segment1>.{9,10})(?P<segment2>CAGAGC){s<=2}(?P<segment3>.{8})(?P<segment4>.{10}).*$
   Read1_src2=^.*(?P<segment1>.{9,10})(?P<segment2>CAGAGC){s<=2}(?P<segment3>.{8})(?P<segment4>.{10}).*$
   Read2_src=(?P<segment5>.*)

In this setting file, 2 read patterns are set for read1. In this case,
first the pattern in Read1\_src are searched and if it is not found,
then the pattern in Read1\_src2 is used for search.

.. _header-n231:

Simple format for read structure specification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since the regular expression matching is computationaly slow, 'simple'
format can also be available for the read structure specification.
Especially for simple structure with fixed component length and no
constant sequence, this 'simple' format will be better than regular
expression. Example for the 10X Chromium v3:

.. code:: shell

   [import]
   src_raw_components=segment1, segment2, segment3
   barcodes=segment1,segment2
   Read1_src=<segment1:16><segment2:12>
   Read2_src=<segment3:0>

Here, simplified read structure is shown in Read1\ *src and Read2*\ src
where 16 bases for segment1 (In this case, cell barcode) and 12 bases
for segment 2 (UMI) are sequenced in read1 and all bases of read2 are
sequenced as segment3 (biological sequence).

**Note:**\ The length '0' means the segment from the current position to
the end. For example,

.. code:: 

   <segment1:10><segment2:15><segment3:0>

means the first 10 bases for segment1, next 15 bases for segment2 and
the lest bases are assigned to segment3.

.. _header-n245:

Computational time and memory consumption
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Computational time linearly depends on the size of input file size.
Memory consumption depends on the sequence length. Basically, it
requires 2~5GB memory for 4million reads processing. For 4 million
reads, sciRNA-seq3(R1:36bp, R2:56bp) barcode extraction with regular
expression took about 9mins and 1.5GB memory, and 10X Visium(R1:28bp,
R2:120bp) barcode extraction with simple format took about 12mins and
3.5GB memory.

When you use barista in cluster server system, splitting input fastq
file into 1million~4million reads per file is recommended.
