.. _header-n9:

Quick start
===========

In this practice, we demonstrate the fastq format conversion of
sciRNA-seq3(\ `Cao et al.,
2019 <https://www.nature.com/articles/s41586-019-0969-x>`__) into 10X
Chromuim V3. Small dataset (10000 reads) of sciRNA-seq3 is stored in
barista/quickstart/dataset. Each step requires ~2min and ~5G memory.

To see the available arguments and the grammer of config file in detail,
see each command page.

-  import<link!>

-  correct<link!>

-  convert<link!>

-  export<link!>

.. _header-n27:

Extract barcodes
----------------

First, fastq file is parsed by ``import``. This ``import`` command
accepts various sequence structure which can be presented in extended
regular expression (fuzzy matching) allowing mismatch.

.. code:: shell

   barista \
   -cmd import \
   -conf ./quickstart/config/example.conf \
   -R1 ./quickstart/dataset/sci3_small_R1.fastq \
   -R2 ./quickstart/dataset/sci3_small_R2.fastq \
   -d ./quickstart/out_import \
   -o barista_import

.. code:: shell

   #output directory structure
   ./quickstart/out_import/
   ├── barista_import_srcCount.pkl.gz
   ├── barista_import_srcQual.tsv.gz
   └── barista_import_srcSeq.tsv.gz

You can see the extracted sequences in Seq.tsv.gz file.

.. _header-n46:

Correct barcodes
----------------

Extracted sequences are corrected by ``correct``. ``correct`` command
adopts 2-step correction, minority correction and whitelist correction
(if whitelist available). Minority correction is based on an assumption
that the true barcodes are more detected than erroneous barcodes.
Therefore, the minority barcode pool are compared with the majority pool
and corrected to the most feasible barcode with a given Levenstein
distance. If the sequencing platform provides pre-defined barcode
sequence list called as 'whitelist', the majority barcodes are compared
with the corresponding whitelist at the 2nd step.

Here, we correct 2 cell barocdes in sciRNA-seq3 with distance=1 for
minority correction and distance=2 for whitelist correction.

.. code:: shell

   barista \
   -cmd correct \
   -conf ./quickstart/config/example.conf \
   -ip ./quickstart/out_import/barista_import_srcCount.pkl.gz \
   -d ./quickstart/out_correct/ \
   -o barista_correct

.. code:: shell

   #output directory structure
   ./quickstart/out_correct/
   ├── barista_correct_barcodeSummary_segment1.png
   ├── barista_correct_barcodeSummary_segment4.png
   ├── barista_correct_srcCorrect.pkl.gz
   └── barista_correct_srcReference.tsv.gz

After correction, barcode rank plot are exported in
./quickstart/out_correct/ directory.

(figure here!)

Then, extracted sequences generated in ``import`` are corrected by
``mk_sval`` command. This ``mk_sval`` command generates corrected
sequence table and barcode value table which will be used in barcode
conversion step.

.. code:: shell

   barista \
   -cmd mk_sval \
   -conf ./quickstart/config/example.conf \
   -crp ./quickstart/out_correct/barista_correct_srcCorrect.pkl.gz \
   -rs ./quickstart/out_import/barista_import_srcSeq.tsv.gz \
   -rq ./quickstart/out_import/barista_import_srcQual.tsv.gz \
   -d ./quickstart/out_mk_sval/ \
   -o barista_mk_sval

.. code:: shell

   #output directory structure
   ./quickstart/out_mk_sval/
   ├── barista_mk_sval_correct_result.tsv.gz
   ├── barista_mk_sval_correct_srcQual.tsv.gz
   └── barista_mk_sval_correct_srcValue.tsv.gz

Extracted sequences after correction are stored in result.tsv file.

.. _header-n216:

Barcode conversion
------------------

For the barcode conversion, 2 types of conversion, 'global convert' and
'local convert', are considered. In global convert, source barcodes are
naively converted into destination barcode. In local convert, barcode
conversion of 1 barcode depends on other barcode(s). For example,
scRNA-seq UMI is meaningful at least under the same cell barcode and UMI
'AAAAA' in cell A and 'AAAAA' in cell B can tag different molecule in
defferent cells.

The first step of conversion is globalConvert:

.. code:: shell

    barista \
    -cmd globalConvert \
    -conf ./quickstart/config/example.conf \
    -sv ./quickstart/out_mk_sval/barista_mk_sval_correct_srcValue.tsv.gz \
    -sq ./quickstart/out_mk_sval/barista_mk_sval_correct_srcQual.tsv.gz \
    -ref ./quickstart/out_correct/barista_correct_srcReference.tsv.gz \
    -d ./quickstart/out_globalConvert/ \
    -o barista_globalConvert

.. code:: shell

   #output directory structure
   ./quickstart/out_globalConvert/
   ├── barista_globalConvert_constraintTree.pkl.gz
   ├── barista_globalConvert_global_destQual.tsv.gz
   └── barista_globalConvert_global_destValue.tsv.gz

cell barcodes 1 and 2 are now globally converted. Then, for the
preparetion of local conversion of UMI, build a correspondence tree to
store the combination of cell barcode and UMI using ``mergeTree``
function.

.. code:: shell

   barista \
   -cmd mergeTree \
   -conf ./quickstart/config/example.conf \
   -lp ./quickstart/out_globalConvert/barista_globalConvert_constraintTree.pkl.gz \
   -d ./quickstart/out_mergeTreeAfter tree merge, 1 pickle file is exported.

.. code:: shell

   #output directory structure
   ./quickstart/out_mergeTree/
   └── barista_mergedTree.pkl.gz

Finally, UMI is converted by ``localConvert`` command:

.. code:: shell

   barista \
   -cmd localConvert \
   -conf ./quickstart/config/example.conf \
   -sv ./quickstart/out_mk_sval/barista_mk_sval_correct_srcValue.tsv.gz \
   -dv ./quickstart/out_globalConvert/barista_globalConvert_global_destValue.tsv.gz \
   -l d_value_UMI \
   -mp ./quickstart/out_mergeTree/barista_mergedTree.pkl.gz \
   -d ./quickstart/out_localConvert/ \
   -o barista_localConvert

.. code:: shell

   #output directory structure
   ./quickstart/out_localConvert/
   └── barista_localConvert_d_value_UMIglobal_destValue.tsv.gz

If your conversion were succeded, you'll get the message
``All conversions were DONE! Proceed export step.``

.. _header-n141:

Read structure reformatting and export
--------------------------------------

For the final step of conversion, use ``export`` command to generate
fastq file with desired sequence structure. Here the converted
sciRNA-seq3 file is exported as 10X v3 format. This step consumes time
and memory compared to previous other steps, but most of them are for
generation of 3million cell barcodes list and 12 nuucleotides UMI list
which is constant even if the input sequence file is much larger than
this practice dataset.

**Note:** Because of the file size limitaion of Github, 10X barcode
whitelist is gzipped. Please decompress the
./quickstart/dataset/3M-february-2018.txt.gz before you run the
following export script.

.. code:: shell

   barista \
   -cmd export \
   -conf ./quickstart/config/example.conf \
   -dv ./quickstart/out_localConvert/barista_localConvert_d_value_UMIglobal_destValue.tsv.gz \
   -dq ./quickstart/out_globalConvert/barista_globalConvert_global_destQual.tsv.gz \
   -rs ./quickstart/out_import/barista_import_srcSeq.tsv.gz \
   -rq ./quickstart/out_import/barista_import_srcQual.tsv.gz \
   -d ./quickstart/out_export/ \
   -o barista_export_10X

.. code:: shell

   #output directory structure
   ./quickstart/out_export/
   ├── barista_export_10XRead1_dest_destSeq.fastq.gz
   └── barista_export_10XRead2_dest_destSeq.fastq.gz

These 10X formatted fastq file can be analyzed by Cellranger.

.. _header-n155:

Read demultiprexing
-------------------

barista can split the input fastq files based on specific barcode using
``demultiplex`` command. ``demultiplex`` can be performed after
generation of result file in ``mk_sval`` step. In this demonstration,
the cell barcode1 is used for demultiplex tag and cell barcode2+UMI and
biological read are exported as Read1 and Read2, respectively.

.. code:: shell

   barista \
   -cmd demultiplex \
   -conf ./quickstart/config/example.conf \
   -cs ./quickstart/out_mk_sval/barista_mk_sval_correct_result.tsv.gz \
   -rq ./quickstart/out_import/barista_import_srcQual.tsv.gz \
   -cq ./quickstart/out_mk_sval/barista_mk_sval_correct_srcQual.tsv.gz \
   -d ./quickstart/out_demultiplex/ \
   -o barista_demulti

.. code:: shell

   #output directory structure
   ./quickstart/out_demultiplex/
   └── barista_demulti_20200414085111
       ├── barista_demulti_AAACCATAGT_20200414085111_Read1.fastq.gz
       ├── barista_demulti_AAACCATAGT_20200414085111_Read2.fastq.gz
       ├── barista_demulti_AAACTCCAAT_20200414085111_Read1.fastq.gz
       ├── barista_demulti_AAACTCCAAT_20200414085111_Read2.fastq.gz
       ├── barista_demulti_AAACTGAACT_20200414085111_Read1.fastq.gz
       ├── barista_demulti_AAACTGAACT_20200414085111_Read2.fastq.gz
       ├── barista_demulti_AAAGCTGATT_20200414085111_Read1.fastq.gz
       ├── barista_demulti_AAAGCTGATT_20200414085111_Read2.fastq.gz
       ...

768 fastq files(384 cell barcode1 \* 2) are generated.

.. _header-n230:

Read tagging
------------

Similar to ``demultiplex`` command, barista also implement ``tag``
command. ``tag`` exports fastq file in which desired barcodes are tagged
on the header line. This command is useful when some barcodes are needed
to be attached to target segment, for example, extract biological read
with cell barcode+UMI in scRNA-seq. Here, we demonstrate to tag cell
barcode 1 and 2 with UMI for read1 and cell barcode 1, 2 and UMI with
constant sequence+biological sequence for read2.

.. code:: shell

   barista 
   -cmd tag 
   -c ./quickstart/config/example.conf 
   -cs ./quickstart/out_mk_sval/barista_mk_sval_correct_result.tsv.gz 
   -rq ./quickstart/out_import/barista_import_srcQual.tsv.gz 
   -cq ./quickstart/out_mk_sval/barista_mk_sval_correct_srcQual.tsv.gz 
   -d ./quickstart/out_tag/ 
   -o barista_tag

.. code:: shell

   #output directory structure
   ./quickstart/out_tag/
   ├── barista_tagRead1.fastq.gz
   └── barista_tagRead2.fastq.gz
