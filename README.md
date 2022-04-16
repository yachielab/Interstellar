# INTERSTELLAR

## Overview

INTERSTELLAR (interpretation, scalable transformation, and emulation of large-scale sequencing reads) is a versatle software tool that extracts data values encoded from sequencing reads and translates them into sequencing reads of another structure according to a user-defied process configuration file.

<center><img src="https://github.com/yachielab/Interstellar/blob/img/image/fig01.jpg" width="750 px"/></center>

It first identifies sequence segments of reads in FASTQ files according to the user’s definition provided in a process configuration file (Step 1). The source read structure can be defined flexibly using regular expression and/or predefined allowlists. A combinatorial segment group can be defined to collectively denote a specific value. A parental segment (or combinatorial parental segment group) can be paired with an independent set of local segments (or combinatorial local segment groups), where sequence-to-value interpretations of local segments are independently defined for each parental segment. For example, cell IDs and UMIs of typical scRNA-seq reads can be defined as parental segments and their local segments, respectively, where the same UMI sequences associated to different cell IDs are interpreted as different objects. Multiple source read structures can also be defined for a single set of input sequencing reads that are produced by a one-shot sequencing of different libraries.

The segment identification process can be performed independently for fragmented FASTQ files using distributed computing, where each fragmented process yields segmented sequences and count-up information for each unique segment sequence. The sequence count-up information derived from different fragmented processes is then merged to compute an error-corrected sequence for each unique segment sequence (Step 2). INTERSTELLAR enables four error correction options: imputation-to-majority, mapping-to-allowlist, Bartender, and a user-developed plugin. Once an error-correction rule table is generated, it is used to error-correct segment sequences originated from each of the fragmented FASTQ files using distributed computing. The above-mentioned read interpretation process can be applied to theoretically any of the high-throughput sequencing read analysis and the generated error-corrected segment sequence files enable efficient development of their downstream data analysis pipelines.

If defined in the process configuration file, the read translation into destination read structures are next processed for the source reads with all segment sequences required for the translation to have identified in the previous processes (Step 3). Destination read structures can flexibly be specified by using IUPAC codes and allowlists of destination segment sequences. First, using distributed computing, a segment value file and a value tree are extracted from each of the error-corrected segment sequence files, where each unique segment sequence is converted into a numerical value and parental-local segment allocations of unique values and unique combinatorial value groups are represented in a tree structure. The value tree files originated from the fragmented FASTQ files are then passed to a single computing node to generate a merged value tree. Next, the values in the merged value tree are replaced by new values to minimize the number of numerical value species for each variable in a way that they still uniquely maintain the same tree topology (Step 4). Obtaining the value conversion rule table that achieved the optimization of the merged value tree, segment value files are separately processed to derive optimized segment value files and then destination FASTQ files by distributed computing, where unique value-to-destination-sequence conversion rules are autonomously generated for the destination read structures (Steps 5 and 6).

Throughout the process, the average Q scores of source segment sequences are bequeathed from the fragmented FASTQ files through the intermediate segment sequence and value files and given to all letters of corresponding destination segments in the generating FASTQ files. New bases that are not associated to the values derived from the source reads are all given a Q score of 40 in the destination reads.

The distributed computing process is designed to perform many small conversion tasks in parallel, where the generation of the conversion rules that require to monitor the entire segment sequence or value space is operated using a single computing node by compressing information from each fragmented task into a small hash table data.



## Installation

### Software dependencies:

- Python3 (version: 3.7.0 or later)
- matplotlib <https://anaconda.org/conda-forge/matplotlib>
- numpy <https://anaconda.org/conda-forge/numpy>
- pandas <https://anaconda.org/conda-forge/pandas>
- regex <https://pypi.org/project/regex/>
- kneed <https://anaconda.org/conda-forge/kneed>
- symspellpy <https://pypi.org/project/symspellpy/> 

The following versions have been tested for INTERSTELLAR 1.0.0:

- Python 3.7.3
- matplotlib 3.4.2
- numpy 1.20.2
- pandas 1.2.4
- regex 2021.4.4
- kneed 0.7.0
- symspellpy 6.7.0

The following tools can optionally be used. Set the paths to $PATH.

- Seqkit <https://bioinf.shenwei.me/seqkit/> (Required for the distributed computing mode)
- FLASH2 <https://github.com/dstreett/FLASH2> (Required for the merging of paired-end reads)
- Bartender <https://github.com/> (Required for the error correction using Bartender) 

### Installation

Obtain INTERSTELLAR by:   	

```bash
git clone /URL/to/Interstellar.git Add the executable authority to Interstellar/py/*
```

 and set the `Interstellar/py/` directory to $PATH:

```bash
cd Interstellar
chmod u+x py/*
```

 Confirm if all the libraries are successfully imported:

```bash
Interstellar -h 
```

### Distributed computing

To enable distributed computing, `qsub` and `qacct` commands must be available to submit and monitor jobs. We have tested the distributed computing mode with UGE 8.6.12 and SGE 6.2u6.



## Execution
### General usages

INTERSTELLAR accepts two different configuration files, one provided with -conf (required) to configure sequencing read interpretation and translation and the other provided with -qconf (optional) to configure distributed computing setting. When a distributed computing is configured, INTERSTELLAR splits input FASTQ files into small fragmented files and processes them separately.

To run INTERSTELLAR with a single computing node (no distributed computing)

```bash
Interstellar -conf process.conf
```

To run INTERSTELLAR with a distributed computing mode:

```bash
Interstellar -conf process.conf -qconf distributed_computing.conf
```

See below for the detail settings of the configuration file.

### Outputs

INTERSTELLAR first performs a process named `value_extraction`, and then performs `value_translation` and/or  `demultiplex` based on the user’s definition in the configuration file. All INTERSTELLAR processes are performed in a user-defined project directory. Output files produced by each process are exported to `out` child directory under the process directory in the main project directory. An example project directory structure follows:

```
./Project_directory/
└── sample1/
    ├── demultiplex/
    │       └── out/
    ├── value_extraction/
    │       └── out/
    └── value_translation/
            └── out/
```

Following files are produced by these processes:

| **Step name**     | **Output files**                                             |
| ----------------- | ------------------------------------------------------------ |
| value_extraction  | A TSV file for identified sequence segments that are quality-filtered and error-corrected . |
| value_translation | FASTQ files for translated reads.                            |
| demultiplex       | FASTQ or TSV files demultiplexed according to user’s definitions. |



## Configuration files

### Process configuration file

A process configuration file consists of process setting sections, `general`, and `value_extraction` to configure how to interpret input reads and `value_translation` or `demultiplex` to configure how to handle the interpreted information. 

- general (required)
  In this section, users configure basic information required for INTERSTELLAR, such as the project directory path.
- value_extraction (required)
  In this section, users configure how sequencing reads are interpreted. This section produces a identified sequence segment table in a TSV format. The process configure in this section is necessary for any of its downstream processes, value_translation, demultiplex, or annotate_header.
- value_translation (optional)
  In this section, users configure value translation process and destination read structure. Source segment values extracted in value_extraction are optimized, considering the value combinations and parental-local segment allocations before the translation into destination segment sequences. The desination segment sequences are concatnated and exporterd into output FASTQ files.
- demultiplex (optional)
  In this section, instead of executing value_translation users configure a process to demultiplex the input reads according to the information extracted at value_extraction and annotate the headers of the demultiplexed reads. 

### Distributed computing configuration file

This file consists of a single setting section "qsub" to configure the settings for distributed computing.

In both of the configuration file types, each setting section must start by [section name] to indicate section blocks as follows:

```
[general]
xxx=...
...
[value_extraction]
xxx=...
...
[value_translation]
xxx=...
...
[demultiplex]
xxx=...
...
```



#### Process configuration: general

**Example**

```
[general]
PROJECT_DIR=/path/to/project_directory
SET_SHELL_ENV=/path/to/set_shell_env.sh
SAMPLESHEET=/path/to/samplesheet.tsv
```

**Description**

- PROJECT_DIR (required)
  Absolute path to a project directory. All processes will be performed under this directory. If INTERSTELLAR is executed without a directory in the designated path, a new directory will be created. When a directory exists in the designated path, INTERSTELLAR seeks for completed processes in the configuration file and resumes the remaining processes if any. For example, if a command error happens at value_translation, INTERSTELLAR skips the value_extraction step when it is re-executed with a modified value_translation configuration.

- SET_SHELL_ENV (required)
  Path to a shell script to activate your environment. Example file (https://github.com/yachielab/Interstellar/blob/main/example-dataset/templates/set_shell_env).

- SAMPLESHEET (required)
  A TSV file to specify the correspondence between the input FASTQ file prefixes and their sample names. The input FASTQ file names need to contain the prefix, like \<fileprefix>_R1_001.fastq (or fastq.gz). This allows INTERSTELLAR to interpret dinsinct samples independently. 

  Example samplesheet format:

  ```
  svf_plate10_SRR5664329		plate10
  svf_plate03_SRR5664331		plate03
  svf_plate04_SRR5664330		plate04
  ```

  For this project, the input FASTQ file names can be like:

  ```bash
  /path/to/Input_dir/R1/
  ├── svf_plate10_SRR5664329_R1-001.fastq.gz #sample group: plate10
  ├── svf_plate10_SRR5664329_R1-002.fastq.gz #sample group: plate10
  ├── svf_plate03_SRR5664331_R1-001.fastq.gz #sample group: plate03
  ├── svf_plate03_SRR5664331_R1-002.fastq.gz #sample group: plate03
  ├── svf_plate04_SRR5664330_R1-001.fastq.gz #sample group: plate04
  └── svf_plate04_SRR5664330_R1-002.fastq.gz #sample group: plate04
  
  /path/to/Input_dir/R2/
  ├── svf_plate10_SRR5664329_R2-001.fastq.gz #sample group: plate10
  ├── svf_plate10_SRR5664329_R2-002.fastq.gz #sample group: plate10
  ├── svf_plate03_SRR5664331_R2-001.fastq.gz #sample group: plate03
  ├── svf_plate03_SRR5664331_R2-002.fastq.gz #sample group: plate03
  ├── svf_plate04_SRR5664330_R2-001.fastq.gz #sample group: plate04
  └── svf_plate04_SRR5664330_R2-002.fastq.gz #sample group: plate04
  ```



#### Process configuration: value_extraction

**Example**

```
[value_extraction]
# Input file path
READ1_PATH=/path/to/Input_dir/R1/*_R1-*.fastq.gz
READ2_PATH=/path/to/Input_dir/R2/*_R2-*.fastq.gz

# Read structure configuration
READ1_STRUCTURE=^(?P<segment1>.{9,10})(?P<segment2>CAGAGC){s<=2}(?P<segment3>.{8})(?P<segment4>.{10}).*$
READ2_STRUCTURE=(?P<segment5>.*)

# Parental-local variable allocations
segment1,segment3=CHILD_OF(segment1,segment4)

# Quality filtering and correction for each sequence segment
cb_value1=QUALITY_FILTER(source:segment1,min_nucleotide_Q-score:5,min_avg_Q-score:20)  >> M2A_CORRECTION(levenshtein_distance:1,path:/path/to/Intersellar/example-dataset/sci-RNA-seq3_HP_purified_allowlist.txt) 
>> SEQ2VALUE()

cb_value2=QUALITY_FILTER(source:segment4, min_nucleotide_Q-score:5, min_avg_Q-score:20) >> I2M_CORRECTION(rank:auto, levenshtein_distance:1) >>  SEQ2VALUE()

umi_value3=QUALITY_FILTER(source:segment3, min_nucleotide_Q-score:5, min_avg_Q-score:20) >>  SEQ2VALUE()

bio_value5=QUALITY_FILTER(source:segment5, min_nucleotide_Q-score:5, min_avg_Q-score:20)
```

**Description**

 The following variable names are fixed for INTERSTELLAR: 

| Name                 | Description                                                  |
| -------------------- | ------------------------------------------------------------ |
| READ1_PATH           | File directory path or file paths for Read 1 files           |
| READ2_PATH           | File directory path or file paths for Read 2 files           |
| INDEX1_PATH          | File directory path or file paths for Index 1 files          |
| INDEX2_PATH          | File directory path or file paths for Index 2 files          |
| FLASH                | Read pairs to be merged by FLASH2. The value must be either "READ1-READ2”, "READ1-INDEX1”, "READ1-INDEX2”, "READ2-INDEX1”, "READ2-INDEX2”, or "INDEX1-INDEX2”. Currently, you can use FALSH only for a single pair of reads. |
| FLASH_MIN_OVERLAP    | Minimum overlap length requied for the read merging by FLASH2. |
| FLASH_MAX_OVERLAP    | Maximum overlap length allowed for the read merging by FLASH2. |
| READ1_STRUCTURE      | Read structure(s) of Read 1 reads defined by regular expression. Multiple structures can be given like "READ1_STRUCTURE_1=...”, "READ1_STRUCTURE_2=...”, … |
| READ2_STRUCTURE      | Read structure(s) of Read 2 reads defined by regular expression. Multiple structures can be given like "READ2_STRUCTURE_1=...”, "READ2_STRUCTURE_2=...”, … |
| INDEX1_STRUCTURE     | Read structure(s) of Index 1 reads defined by regular expression. Multiple structures can be given like "INDEX1_STRUCTURE_1=...”, "INDEX1_STRUCTURE_2=...”, … |
| INDEX2_STRUCTURE     | Read structure(s) of Index 2 reads defined by regular expression. Multiple structures can be given like "INDEX2_STRUCTURE_1=...”, "INDEX2_STRUCTURE_2=...”, … |
| READ_FLASH_STRUCTURE | Read structure of merged reads by FLASH2 defined by regular expression. Multiple structures can be given like "READ_FLASH_STRUCTURE_1=...”, "READ_FLASH_STRUCTURE_2=...”, … |

Variable names of identified sequence segments can be freely defined by the user.

- READ1_PATH, READ2_PATH, INDEX1_PATH, INDEX2_PATH

  Path to input FASTQ files. FASTQ file names must include the fleprefixes indicated in SAMPLESHEET.
  Example:

  ```
  READ1_PATH=/path/to/Input_dir/R1/*.fastq.gz
  READ2_PATH=/path/to/Input_dir/R2/*.fastq.gz
  ```

  Alternatively, user can provide a directory path containing FASTQ files.

  ```
  READ1_PATH=/path/to/Input_dir/R1/
  READ2_PATH=/path/to/Input_dir/R2/
  ```

- FLASH

  Paired-end reads that sequenced overlapping regions can be merged into single reads. Only a single pair of two reads among READ1, READ2, INDEX1 and INDEX2 can be merged. The pair should be specified as "READ1-READ2”, "READ2-READ1”, "READ1-INDEX1”, "INDEX1-READ1”, "READ1-INDEX2”, "INDEX2-READ1”, "READ2-INDEX1”, "INDEX1-READ2”, "READ2-INDEX2”, "INDEX2-READ2”, "INDEX1-INDEX2"or "INDEX2-INDEX2”. Direction of the merged read will follow that of the first read to be merged.

  Example configurations for FLASH are shown in Example 5(https://github.com/yachielab/Interstellar/blob/main/example-dataset/Example5/ex5.conf).

- FLASH_MIN_OVERLAP
  Corresponding to the `-m` option of FLASH2 thatspecifies the minimum overlap length. Default=20.

- FLASH_MAX_OVERLAP
  Corresponding to the `-M` option of FLASH2 that specifies the maximum overlap length. Default=30.

- READ1_STRUCTURE, READ2_STRUCTURE, INDEX1_STRUCTURE, INDEX2_STRUCTURE, READ_FLASH_STRUCTURE, …

  Read structures specified by regular expression. Multiple read structures can be specified for a single read, like READ1_STRUCTURE_1=..., READ1_STRUCTURE_2=..., … or READ_FLASH_STRUCTURE_1=..., READ_FLASH_STRUCTURE_2=..., … 

  Variables specified by ?P<> are user-defined sequence segments. INTERSTELLAR adopts the extended regular expression implemented in the regex library (https://pypi.org/project/regex/) that allows mismatches.

  

  Ex.1) Quartz-seq library

  ```
  READ1_STRUCTURE=^(?P<seg1>.{14})(?P<seg2>.{8})$
  READ2_STRUCTURE=^(?P<seg3>.*)$
  ```

  In this example, READ1 is 22-bp long in total and consists of a 14-bp "seg1"segment and an 8-bp "seg3"segment. The entire sequence of READ2 is obtained as a "seg3"segment. Here, seg1, seg2, and seg3 are a cell ID, UMI, and cDNA, respectively.

  

  Ex.2) sci-RNA-seq3 library

  ```
  READ1_STRUCTURE=^(?P<seg1>.{9,10})(?P<seg2>CAGAGC){s<=2}(?P<seg3>.{8})(?P<seg4>.{10})$
  ```

  Here we demonstrate how to express Read 1 of a sci-RNA-seq3 library. In a 34- or 33-bp sequence, a constant sequence of "CAGAGC"(seg2) is allowed to have up to two mismatches and sandwithced by 9- or 10-bp "seg1"on the left and 8-bp "seg3"and 10-bp "seg4"on the right. Here, seg1, seg2, seg3, and seg4 are a cell ID conferred by ligation, constant sequence, UMI, and another cell ID conferred by reverse transcription.

- Parental-local value allocations

  Parental-loval value allocations of the segments can be defined by `CHILD_OF()`. 

  Ex. 1)  seg1 is a local variable of seg2:

  ```
  seg1 = CHILD_OF(seg2)
  ```

  Ex. 2)  seg2 and seg3 together consist of a combinatorial variable where seg1 is local to them:

  ```
  seg1 = CHILD_OF(seg2,seg3)
  ```

  Ex. 3) a combinatorial variable (seg1, seg2) is local to another combinatorial variable (seg3, seg4):

  ```
  seg1,seg2 = CHILD_OF(seg3,seg4)
  ```

- Quality filtration, sequence error-correction (four different options), and sequence-to-value conversion

  - QUALITY_FILTER()

    Filter segment sequences according to their Q scores with the following arguments:

    - source (required)

      Sequence segment variable name.

    - min_nucleotide_Q-score

      Minimum per-base Q score threshold. Default=5.

    - min_avg_Q-score

      Minimum average Q score threshold.. Default=20.

  - I2M_CORRECTION()

    Error-correct segment sequences by "imputation-to-majority."

    - source

      Sequence segment varialble name. Not required If defined upstream and inherirted by ">>"(see below).

    - rank

      Rank threshold for the minority-majority boundary of the rank-count plot of unique segment sequences. If this set to "auto," the knee-point of the rank-count plot is adopted. Default=auto.

    - levenshtein_distance

      Maximum Levenshtein distance threshold allowed to correct sequences in the minority group to their nearest sequences in the majority group. Default=1.

  - M2A_CORRECTION()

    Error-correct segment sequences by “mapping-to-allowlist.”

    - source

      Sequence segment variable name. Not required if defined upstream and inherited by “>>” (see below).

    - path (required)

      Path to the sequence allowlist.

    - levenshtein_distance

      Maximum Levenshtein distance threshold allowed to correct the input sequences to their nearest sequences in the allowlist. Default=1.

  - BARTENDER_CORRECTION()

    Error-correct segment sequences using Bartender (URL). Please see Bartender’s manual for details. 

    - source

      Sequence segment variable name. Not required if defined upstream and inherited by “>>” (see below).

    - -c

      Sequence frequency cutoff.

    - -d

      Maximum cluster distance that may be merged.

    - -z

      Cluster merging threshold.

    - -l

      Seed length.

    - -s

      Number of non-overlapping positions between two adjacent seeds.

  - CUSTOM_CORRECTION()

    Error-correct segment sequences using a user-provided shell script.

    - source

      Sequence segment variable name. Not required if defined upstream and inherited by “>>” (see below).

    - shell_script

      Path to the shell script.

    The script must be designed to accept an input sequence segment list in a TSV (tab-separated value) file format as the first argument $1 and output the conversion result in STDOUT in a TSV format at the second argument $2 of the shell script.

    

    Input data file format:

    ```
    AAGGC	0
    GAGGT	1
    GCGTA	2
    AAGGC	3
    …
    ```

    The first and second column stores extracted segment sequences and the line numbers are shown in the second column. Each line represents each FASTQ record after quality filtering and the sequences are not de-duplicated. 

    

    Output data format:

    ```
    AAGGC	AAGGC
    GAGGT	GAGGG
    GCGTA	-
    AAGGC	AAGGC
    ```

    See the example shell script.

  - SEQ2VALUE()

    Convert segment sequences into (unoptimized) segment values. This process is required if the segment sequences are supposed to be translated into different sequences.

    - source
      Sequence segment variable name. Not required if defined upstream and inherited by “>>” (see below).

  - Pipe-like operator

    “>>” can be used as like a pipe like operator, connecting multiple processes in one line.

    

  Ex. 1) Quality filteration of “segment1” followed by the mapping-to-allowlist correction and the sequence-to-value conversion.:

  ```
  segment1.filtered.corrected.value=QUALITY_FILTER(source:segment1,min_nucleotide_Q-score:5,min_avg_Q-score:20) 
  >> M2A_CORRECTION(levenshtein_distance:1,path:./path/to/allowlist.txt) 
  >>  SEQ2VALUE()
  ```

  Ex. 2) The same operation can be described as follows without the “>>” operator.:

  ```
  segment1.filtered=QUALITY_FILTER(source:segment1,min_nucleotide_Q-score:5,min_avg_Q-score:20) segment1.filtered.corrected=M2A_CORRECTION(source: segment1.filtered, levenshtein_distance:1,path:/path/to/allowlist.txt) segment1.filtered.corrected.value=SEQ2VALUE(source:segment1.filtered.corrected)
  ```



#### Process configuration: value_translation

**Example**

```
[value_translation]
#Value to destination sequence conversion
dest_segment1=VALUE2SEQ(source:cb_value1,cb_value2, allowlist_path:/path/to/allowlist.txt)
dest_segment2=VALUE2SEQ(source:umi_value3, randseq_pattern:NNNNNNNNNNNN)
dest_segment3=bio_value5
exConst1=”TCTC”

#Read structure configuration
READ1_STRUCTURE=dest_segment1+dest_segment2+exConst1
READ2_STRUCTURE=dest_segment3
```

**Desctiprion**

The following variable names are fixed for INTERSTELLAR:

| Name             | Description                                                  |
| ---------------- | ------------------------------------------------------------ |
| READ1_STRUCTURE  | Destination read structure of READ1 defined by a sequence of destination segments connected using “+”. |
| READ2_STRUCTURE  | Destination read structure of READ2 defined by a sequence of destination segments connected using “+”. |
| INDEX1_STRUCTURE | Destination read structure of INDEX1 defined by a sequence of destination segments connected using “+”. |
| INDEX2_STRUCTURE | Destination read structure of INDEX2 defined by a sequence of destination segments connected using “+”. |

Variable names of destination sequence segments can be freely defined by the user.

- Defining destination segment sequences

  The destination sequence segments can be defined by VALUE2SEQ() or SEQ2SEQ(), or by directly assigning source sequence segments. A single destination segment can be defined as follows:

  ```
  dest = VALUE2SEQ()
  ```

  Combinatorial destination segments corresponding to a single source value (like in the cell ID conversion from 10X Chromium to SPLiT-seq) can be defined as follows: 

  ```
  dest1,dest2 = VALUE2SEQ(...)
  ```

  The following functions can be used for producing defining sequence segments.

  - VALUE2SEQ()

    Assign destination segment sequences to values so they fit to an allowlist or destination sequence segment structure. In this section, value space optimization is automatically performed for source segment values according to the segment grouping rules and parental-local segment value allocations defined in the value_extraction section. The value space optimization is not performed if any parental-local segment value allocation is defined earlier.

    - source (required)

      Value names. Comma-separated values if multiple vlaules are translated into a sequence segment. Note that the source segment combinations used in 

    - allowlist_path

      Path to the destination sequence allowlist. If multiple destination sequence segments are to be generated, provide their corresponding allowlists as comma-separated values in the same order of destination sequence segments.

    - randseq_pattern

      Destination sequence segment structure defined by IUPAC codes. If multiple destination sequence segments are to be generated defined, provide their corresponding allowlists as comma-separated values in the same order of destination sequence segments. 

  - SEQ2SEQ()
    Assigns destination segment sequences defined for their corresponding source segment sequences according to a user-defined sequence conversion table.

    - source (required):

      Source sequence segment name. Comma-separated values if multiple source segments are converted into a destianation sequence(s).

    - conversion_table

      Path to the sequence conversion table (TSV file format). The header line must correspond to the target source and destination sequence segment names defined in the configuration file. Each SEQ2SEQ function requires a unique sequence conversion table.

      Ex. 1: A sequence conversion table to convert specific combinations of source1 and source2 into their corresponding combinations of dest1 and dest2.

      ```
      source1		  source2		  dest1		      dest2
      ACTTGATTGT	TCCTACCAGT	ACTGAGTGC    	CGAGATGT
      ACTTGATTGT	GCGTTGGAGC	TGATCGTCATG  	GCGCATTC
      ACTTGATTGT	GATCTTACGC	GATCGGTACG   	AACCCTTG
      ACTTGATTGT	CTGATGGTCA	GGCCCTTA     	CGGACAAC
      ACTTGATTGT	CCGAGAATCC	ATCTGCATC    	AGTTTAGA
      ACTTGATTGT	GCCGCAACGA	GATACAAACT  	ATATACCT
      ACTTGATTGT	TGAGTCTGGC	TTATAGCC     	GGTGGGAT
      ACTTGATTGT	TGCGGACCTA	AATTCAGTA    	AAGGTAAT
      ACTTGATTGT	ACGGAGGCGG  GATCCGTCCA   	ATGTTGGC
      ```

      

  - Bequeathing a source segment
    Source sequence segment variables can be used to define destination reads.

- READ1_STRUCTURE, READ2_STRUCTURE, INDEX1_STRUCTURE, and INDEX2_STRUCTURE 

  Destination read structures can be defined by concatnating destination sequence segments using “+”.



#### Process configuration: demultiplex

### Example

```
[demultiplex]
FORMAT=FASTQ
KEY=seg1,seg2
READ1_HEADER_ADDITION = seg1,seg2
READ1_STRUCTURE       = seg1+seg3
READ2_HEADER_ADDITION = seg1,seg2
READ2_STRUCTURE       = seg4+seg5
```

**Description**

This section enables demultiplexing of input reads according to sequence segments identified in the [value_extraction] section. 

- FORMAT (required)

  Export data format. FASTQ or TSV.

- KEY

  Source segment(s) to demultiplex the input reads. Use comma-separated values if multiple segments are used for this process. If not provided, the reads are not demultiplexed into separate FASTQ files (only available if FORMAT=FASTQ).

- TARGET

  Source segments to be included in the demultiplexed sequence table. Use comma-separated values for multiple segments. This option is available if FORMAT=TSV.

- READ1_STRUCTURE, READ2_STRUCTURE, INDEX1_STRUCTURE, and INDEX2_STRUCTURE

  Destination read structures can be defined by concatnating destination sequence segments using “+”. This option is available if FORMAT=FASTQ.

- READ1_HEADER_ADDITION, READ2_HEADER_ADDITION, INDEX1_HEADER_ADDITION, and INDEX2_HEADER_ADDITION

  The provided sequence segment(s) will be tagged to the read headers. Multiple sequence segments can be provided by comma separated values. The original header followed by the provided sequence segments are concatnated using “_”. This option is available if FORMAT=FASTQ.



#### Distributed computing configuration: qsub

### Example

```
[qsub]
QOPTION=”-l s_vmem=<int>G,mem_req=<int>G”
MEM_MAX=128
MEM_MIN=6
NUM_READS=2000000
```

**Description**

The distributed computing mode can be configured in this section by defining the following parameteres.

- QOPTION (required)

  Options to be added to the qsub command. This is strongly recommended to be defined as some of the INTERSTELLAR operations require large memory size. In each step of INTERSTELLAR operations, \<int> will be replaced by an empirically estimated sufficient memory size between MEM_MIN and MEM_MAX defined below. The example QOPTION above represents the one can be used in a UGE environment.

- MEM_MAX

  Maximum memory size. Integer or float. The scaling factors (i.e., K, M, G, etc.) should be omitted but provided in QOPTIONS instead as in the above example. Default=128.

- MEM_MIN

  Minimum memory size. Integer or float. The scaling factors (i.e., K, M, G, etc.) should be omitted but provided in QOPTIONS instead as in the above example. Default=6.

- NUM_READS

  A number of reads per split FASTQ file. Default=2000000.



## Example process configuration files

Example 1. Translation of Drop-seq reads for 10X Cell Ranger

Example 2. Translation of sci-RNA3-seq reads for 10X Cell Ranger

Example 3. Translation of Quartz-seq reads into the SPLiT-seq read structure 

Example 4. Translation of Slide-seq reads for 10X Space Ranger

Example 5. Read segmentation and demultiplexing of RCP-PCR reads

All of the example configuration files can be tested on Google Colaboratory. (https://colab.research.google.com/drive/1nuqPK_zQSXFXHu-9gZR5w9EfsQhH6Itl?usp=sharing)



### Example 1. Translation of Drop-seq reads for 10X Cell Ranger

The example process configuration file (Interstellar/example-dataset/Example1/ex1.conf)

```
[general]
# Working directory
PROJECT_DIR=/path/to/project_dir
 
# Path to template shellscript that
SET_SHELL_ENV=/path/to/set_shell_env

SAMPLESHEET=/path/to/example-dataset/Example1/samplesheet.tsv

[value_extraction]
# Input file path
READ1_PATH=/path/to/Interstellar/example-dataset/Example1/R1
READ2_PATH=/path/to/Interstellar/example-dataset/Example1/R2
 
# Read structure patterns by regular expression
READ1_STRUCTURE=^(?P<src_segment1>.{12})(?P<src_segment2>.{8}).*$
READ2_STRUCTURE=^(?P<src_segment3>.*)$
 
# Parental-local variable allocations
src_segment2=CHILD_OF(src_segment1)
 
# Sequence quality filtration
segment1.filtered = QUALITY_FILTER(source: src_segment1, min_nucleotide_Q-score:5, min_avg_Q-score:20)
segment2.filtered = QUALITY_FILTER(source: src_segment2, min_nucleotide_Q-score:5, min_avg_Q-score:20)
segment3.filtered = QUALITY_FILTER(source: src_segment3, min_nucleotide_Q-score:5, min_avg_Q-score:20)
# Knee correction for cell barcodes
segment1.filtered.kneed = I2M_CORRECTION(source:segment1.filtered, levenshtein_distance:1)
# Sequence to value conversion for cell barcodes and UMIs
segment1.value = SEQ2VALUE(source:segment1.filtered.kneed)
segment2.value = SEQ2VALUE(source:segment2.filtered)
 
[value_translation]
# Definition of destination segments
dest_segment1=VALUE2SEQ(source: segment1.value, allowlist_path:/path/to/3M-february-2018.txt)
dest_segment2=VALUE2SEQ(source: segment2.value, randseq_pattern:NNNNNNNNNNNN)
dest_segment3=segment3.filtered
 
# Read structure configuration
READ1_STRUCTURE=dest_segment1+dest_segment2
READ2_STRUCTURE=dest_segment3
```

Execution

```
Interstellar -conf /path/to/ex1.conf
```

Alternatively, the same configuration can be described using the pipe-like operators “>>” as follows:

```
[general]
# Working directory
PROJECT_DIR=/path/to/project_dir

# Path to template shellscript that 
SET_SHELL_ENV=/path/to/set_shell_env

SAMPLESHEET=/path/to/example-dataset/Example1/samplesheet.tsv

[value_extraction]
# Input file path
READ1_PATH=/path/to/Interstellar/example-dataset/Example1/R1
READ2_PATH=/path/to/Interstellar/example-dataset/Example1/R2

# Read structure patterns by regular expression
READ1_STRUCTURE=^(?P<src_segment1>.{12})(?P<src_segment2>.{8}).*$
READ2_STRUCTURE=^(?P<src_segment3>.*)$

# Parental-local value allocations
src_segment2=CHILD_OF(src_segment1)

# Sequence quality filtration
segment1.filtered.kneed.value=QUALITY_FILTER(source: src_segment1, min_nucleotide_Q-score:5, min_avg_Q-score:20) >> I2M_CORRECTION(levenshtein_distance:1) >> SEQ2VALUE()
segment2.filtered.value=QUALITY_FILTER(source: src_segment2, min_nucleotide_Q-score:5, min_avg_Q-score:20) >> SEQ2VALUE()
segment3.filtered =QUALITY_FILTER(source: src_segment3, min_nucleotide_Q-score:5, min_avg_Q-score:20)

[value_translation]
# Definition of destination segments
dest_segment1 = VALUE2SEQ(source: segment1.filtered.kneed.value, allowlist_path:/path/to/3M-february-2018.txt)
dest_segment2 = VALUE2SEQ(source: segment2.filtered.value, randseq_pattern:NNNNNNNNNNNN)
dest_segment3 = segment3.filtered

# Read structure configuration
READ1_STRUCTURE=dest_segment1+dest_segment2
READ2_STRUCTURE=dest_segment3
```



### Example 2. Translation of sci-RNA3-seq reads for 10X Cell Ranger
The example process configuration file (Interstellar/example-dataset/Example2/ex2.conf)

```
[general]
# Working directory
PROJECT_DIR=/path/to/project_dir

# Path to template shellscript
SET_SHELL_ENV=/path/to/set_shell_env

SAMPLESHEET=/path/to/example-dataset/Example2/samplesheet.tsv

[value_extraction]
# Input file path
READ1_PATH=/path/to/Interstellar/example-dataset/Example2/R1
READ2_PATH=/path/to/Interstellar/example-dataset/Example2/R2

# Read structure patterns by regular expression
READ1_STRUCTURE=^(?P<src_segment1>.{9,10})(?P<src_segment2>CAGAGC){s<=2}(?P<src_segment3>.{8})(?P<src_segment4>.{10}).*$
READ2_STRUCTURE=^(?P<src_segment5>.*)$

# Parental-local value allocations
src_segment3=CHILD_OF(src_segment1,src_segment4)

# Sequence filtering and correction for each segment
segment1.filtered.corrected.value =QUALITY_FILTER(source:src_segment1, min_nucleotide_Q-score:5, min_avg_Q-score:20) >> M2A_CORRECTION(levenshtein_distance:1,path:/path/to/Interstellar/example-dataset/allowlists/sci-RNA-seq3_HP_purified_allowlist.txt) >>  SEQ2VALUE()
segment4.filtered.corrected.value =QUALITY_FILTER(source:src_segment4, min_nucleotide_Q-score:5, min_avg_Q-score:20) >> M2A_CORRECTION(levenshtein_distance:1,path:/path/to/Interstellar/example-dataset/allowlists/sci-RNA-seq3_RT_purified_allowlist.txt) >>  SEQ2VALUE()
segment3.filtered.value =QUALITY_FILTER(source:src_segment3, min_nucleotide_Q-score:5, min_avg_Q-score:20) >>  SEQ2VALUE()
segment5.filtered=QUALITY_FILTER(source:src_segment5, min_nucleotide_Q-score:5, min_avg_Q-score:20)

[value_translation]
#Value to destination sequence conversion
dest_segment1=VALUE2SEQ(source:segment1.filtered.corrected.value, segment4.filtered.corrected.value, allowlist_path:/path/to/3M-february-2018.txt)
dest_segment2=VALUE2SEQ(source:segment3.filtered.value, randseq_pattern:NNNNNNNNNNNN)
dest_segment3=segment5.filtered

#Read structure configuration
READ1_STRUCTURE=dest_segment1+dest_segment2
READ2_STRUCTURE=dest_segment3
```

Execution

```
Interstellar -conf /path/to/ex2.conf
```



### Example 3. Translation of Quartz-seq reads into the SPLiT-seq read structure 

The example process configuration file (Interstellar/example-dataset/Example3/ex3.conf)

```
[general]
# Working directory
PROJECT_DIR=/path/to/project_dir

# Path to template shellscript
SET_SHELL_ENV =/path/to/template.sh

SAMPLESHEET=/path/to/example-dataset/Example3/samplesheet.tsv

[value_extraction]
# Input file path
READ1_PATH=/path/to/Interstellar/example-dataset/Example3/R1
READ2_PATH=/path/to/Interstellar/example-dataset/Example3/R2

# Read structure patterns by regular expression
READ1_STRUCTURE=^(?P<src_segment1>.{14})(?P<src_segment2>.{8}).*$
READ2_STRUCTURE=^(?P<src_segment3>.*)$

# Parental-local value allocations
src_segment2=CHILD_OF(src_segment1)

# Sequence filtering and correction for each segment
segment1.filtered.corrected.value =QUALITY_FILTER(source:segment1, min_nucleotide_Q-score:5, min_avg_Q-score:20) >> M2A_CORRECTION(levenshtein_distance:1,path:/path/to/Interstellar/example-dataset/allowlists/CellBarcode_QuartzSeq2_v31.txt) >>  SEQ2VALUE()
segment2.filtered.value=QUALITY_FILTER(source:src_segment2, min_nucleotide_Q-score:5, min_avg_Q-score:20) >>  SEQ2VALUE()
segment3.filtered=QUALITY_FILTER(source:src_segment3, min_nucleotide_Q-score:5, min_avg_Q-score:20)

[value_translation]
#Value to destination sequence conversion
dest_segment1,dest_segment2,dest_segment3=VALUE2SEQ(source: segment1.filtered.corrected.value, allowlist_path:/path/to/Interstellar/example-dataset/allowlists/split.round1.txt,/path/to/Interstellar/example-dataset/allowlists/split.round2.txt,/path/to/Interstellar/example-dataset/allowlists/split.round3.txt)
dest_segment4=VALUE2SEQ(source: segment2.filtered.value, randseq_pattern:NNNNNNNNNN)
dest_segment5=segment3.filtered

#Read structure configuration
READ1_STRUCTURE=dest_segment4+dest_segment3+"GTGGCCGATGTTTCGCATCGGCGTACGACT"+dest_segment2+"ATCCACGTGCTTGAGAGGCCAGAGCATTCG"+dest_segment1
READ2_STRUCTURE=dest_segment5
```

Execution

```
Interstellar -conf /path/to/ex3.conf
```



### Example 4. Translation of Slide-seq reads for 10X Space Ranger
The example process configuration file. (Interstellar/example-dataset/Example4/ex4.conf)

```
[general]
# Working directory
PROJECT_DIR=/path/to/project_dir

# Path to template shellscript
SET_SHELL_ENV=/path/to/template.sh

SAMPLESHEET=/path/to/example-dataset/Example4/samplesheet.tsv

[value_extraction]
# Input file path
READ1_PATH=/path/to/Interstellar/example-dataset/Example4/R1
READ2_PATH=/path/to/Interstellar/example-dataset/Example4/R2

# Read structure patterns by regular expression
READ1_STRUCTURE=^(?P<src_segment1>.{13})(?P<src_segment2>.{8}).*$
READ2_STRUCTURE=^(?P<src_segment3>.*)$

# Sequence filtering and correction for each segment
segment1.filtered.corrected =QUALITY_FILTER(source:src_segment1, min_nucleotide_Q-score:5, min_avg_Q-score:20) >> M2A_CORRECTION(levenshtein_distance:1,path:/path/to/Interstellar/example-dataset/allowlists/slide_cereb1_bc.txt)
segment2.filtered =QUALITY_FILTER(source:src_segment2, min_nucleotide_Q-score:5, min_avg_Q-score:20)
segment3.filtered=QUALITY_FILTER(source:src_segment3, min_nucleotide_Q-score:5, min_avg_Q-score:20)

[value_translation]
#Value to destination sequence conversion
dest_segment1=SEQ2SEQ(source:segment1.filtered.corrected,conversion_table:/path/to/Interstellar/example-dataset/allowlists/slide_to_visium_cereb1.correspo.txt)
dest_segment2=segment2.filtered
dest_segment3=segment3.filtered

#Read structure configuration
READ1_STRUCTURE=dest_segment1+dest_segment2+"AAAA"
READ2_STRUCTURE=dest_segment3
```

Execution

```
Interstellar -conf /path/to/ex4.conf
```



### Example 5. Read segmentation and demultiplexing of RCP-PCR reads (FORMAT=TSV)
The example process configuration file. (Interstellar/example-dataset/Example5/ex5.conf)

```
[general]
# Working directory
PROJECT_DIR=/path/to/project_dir

# Path to template shellscript
SET_SHELL_ENV=/path/to/template.sh

SAMPLESHEET=/path/to/example-dataset/Example5/samplesheet.tsv

[value_extraction]
# Input file path
READ1_PATH=/path/to/Interstellar/example-dataset/Example5/R1
READ2_PATH=/path/to/Interstellar/example-dataset/Example5/R2

# Paired-end reads to be merged by FLASH2. Skip if you don't need to merge.
FLASH=READ1-READ2

# Minimum overlap length for the FLASH2 option "-m". Skip if you don't need to merge.
FLASH_MIN_OVERLAP=20

# Maximum overlap length for the FLASH2 option "-M". Skip if you don't need to merge.
FLASH_MAX_OVERLAP=100

# If reads are merged by FLASH2 and searching for three different read structures
READ_FLASH_STRUCTURE_1=^.+(?P<plate1>.{9})(?P<PS1>TAACTTACGGAGTCGCTCTACG){s<=1}(?P<row_idx>.{9})(?P<DBU1>CCATACGAGCACATTACGGG){s<=2}(?P<bc1>.+)(?P<DBU2>CTAACTCGCATACCTCTGATAAC){s<=2}(?P<lox2272>.+)(?P<DBlox2272>TGTCAGCACTCTGTCAAAATAGATCG){s<=2}....(?P<DBD1>TCGATAGGTGCGTGTGAAGG){s<=2}(?P<bc2>.+)(?P<DBD2>CCTCAGTCGCTCAGTCAAG){s<=2}(?P<col_idx>.{9})(?P<PS2>CAGGACCTAAAGAATCCCATCC){s<=1}(?P<plate2>.{9}).*$
READ_FLASH_STRUCTURE_2=^.+(?P<plate1>.{9})(?P<PS1>TAACTTACGGAGTCGCTCTACG){s<=1}(?P<row_idx>.{9})(?P<DBloxP>TTAAGTAAGTAATGGCAAACAGCTATTAT){s<=2}(?P<loxp>.+){s<=2}(?P<DBU1>CCATACGAGCACATTACGGG){s<=2}(?P<bc1>.+)(?P<DBU2>CTAACTCGCATACCTCTGATAAC){s<=2}(?P<lox2272>.+){s<=2}(?P<DBlox2272>TGTCAGCACTCTGTCAAAATAGATCG){s<=2}(?P<col_idx>.{9})(?P<PS2>CAGGACCTAAAGAATCCCATCC){s<=1}(?P<plate2>.{9}).*$
READ_FLASH_STRUCTURE_3=^.+(?P<plate1>.{9})(?P<PS1>TAACTTACGGAGTCGCTCTACG){s<=1}(?P<row_idx>.{9})(?P<ADU1>CCCTTAGAACCGAGAGTGTG){s<=2}(?P<bc1>.+)(?P<ADU2>CACTCCGTTCGTCACTCAATAA){s<=2}...(?P<ADloxP>CAGCACTCTGTCAAAATAGATCGGAAA){s<=2}(?P<loxp>.+){s<=2}(?P<ADD1>CTCCAGGGTTAGGCAGATG){s<=2}(?P<bc2>.+)(?P<ADD2>CAATCGCACTATCCCGCTG){s<=2}(?P<col_idx>.{9})(?P<PS2>CAGGACCTAAAGAATCCCATCC){s<=1}(?P<plate2>.{9}).*$
READ_FLASH_STRUCTURE_4=^.+(?P<plate1>.{9})(?P<PS1>TAACTTACGGAGTCGCTCTACG){s<=1}(?P<row_idx>.{9})(?P<ADloxP>CAGCACTCTGTCAAAATAGATCGG){s<=2}...(?P<loxp>.+){s<=2}(?P<ADD1>CTCCAGGGTTAGGCAGATG){s<=2}(?P<bc2>.+)(?P<ADD2>CAATCGCACTATCCCGCTG){s<=2}(?P<lox2272>.+){s<=2}(?P<ADlox2272>TTATGGGTATTATGGGTTAGGGATAACAGGGTAATAAGCTCTAAGTAAGTAACGGCCGCC){s<=3}(?P<col_idx>.{9})(?P<PS2>CAGGACCTAAAGAATCCCATCC){s<=1}(?P<plate2>.{9}).*$

# Sequence filtering and correction for each segment
plate1.filtered.corrected=QUALITY_FILTER(source:plate1,min_nucleotide_Q-score:10,min_avg_Q-score:20) >> M2A_CORRECTION(levenshtein_distance:1, path:/path/to/example-dataset/allowlists/rcppcr_illumina_idx.txt)
plate2.filtered.corrected =QUALITY_FILTER(source:plate2,min_nucleotide_Q-score:10,min_avg_Q-score:20) >> M2A_CORRECTION(levenshtein_distance:1, path:/path/to/example-dataset/allowlists/rcppcr_illumina_idx.txt)
row_idx.filtered.corrected =QUALITY_FILTER(source:row_idx,min_nucleotide_Q-score:10,min_avg_Q-score:20) >> M2A_CORRECTION(levenshtein_distance:1, path:/path/to/example-dataset/allowlists/rcppcr_illumina_idx.txt)
col_idx.filtered.corrected =QUALITY_FILTER(source:col_idx,min_nucleotide_Q-score:10,min_avg_Q-score:20) >> M2A_CORRECTION(levenshtein_distance:1, path:/path/to/example-dataset/allowlists/rcppcr_illumina_idx.txt)

bc1.filtered=QUALITY_FILTER(source:bc1,min_nucleotide_Q-score:10,min_avg_Q-score:20)
bc2.filtered=QUALITY_FILTER(source:bc2,min_nucleotide_Q-score:10,min_avg_Q-score:20)
lox2272.filtered=QUALITY_FILTER(source:lox2272,min_nucleotide_Q-score:10,min_avg_Q-score:20)
DBlox2272.filtered=QUALITY_FILTER(source:DBlox2272,min_nucleotide_Q-score:10,min_avg_Q-score:20)
DBloxP.filtered=QUALITY_FILTER(source:DBloxP,min_nucleotide_Q-score:10,min_avg_Q-score:20)
ADloxP.filtered=QUALITY_FILTER(source:ADloxP,min_nucleotide_Q-score:10,min_avg_Q-score:20)
loxp.filtered=QUALITY_FILTER(source:loxp,min_nucleotide_Q-score:10,min_avg_Q-score:20)
ADlox2272.filtered=QUALITY_FILTER(source:ADlox2272,min_nucleotide_Q-score:10,min_avg_Q-score:20)

[demultiplex]
# Clean segment names separating the file.
KEY=plate1.filtered.corrected,plate2.filtered.corrected,row_idx.filtered.corrected,col_idx.filtered.corrected

# Export data format, FASTQ or TSV
FORMAT=TSV

# Clean setgment names to be included in the demultiplexed sequence table if FORMAT=tsv.
TARGET=bc1.filtered,bc2.filtered,lox2272.filtered,DBlox2272.filtered,DBloxP.filtered,ADloxP.filtered,loxp.filtered,ADlox2272.filtered
```

Execution

```
Interstellar -conf /path/to/ex5.conf
```



### Example 6. Read segmentation and demultiplexing of RCP-PCR reads (FORMAT=FASTQ)

The example process configuration file. (Interstellar/example-dataset/Example5/ex5.conf)

```
[general]
# Working directory
PROJECT_DIR=/path/to/project_dir

# Path to template shellscript
SET_SHELL_ENV=/path/to/template.sh

SAMPLESHEET=/path/to/example-dataset/Example6/samplesheet.tsv

[value_extraction]
# Input file path
READ1_PATH=/path/to/Interstellar/example-dataset/Example5/R1
READ2_PATH=/path/to/Interstellar/example-dataset/Example5/R2

# Paired-end reads to be merged by FLASH2. Skip if you don't need to merge.
FLASH=READ1-READ2

# Minimum overlap length for the FLASH2 option "-m". Skip if you don't need to merge.
FLASH_MIN_OVERLAP=20

# Maximum overlap length for the FLASH2 option "-M". Skip if you don't need to merge.
FLASH_MAX_OVERLAP=100

# If reads are merged by FLASH2 and searching for three different read structures
READ_FLASH_STRUCTURE_1=^.+(?P<plate1>.{9})(?P<PS1>TAACTTACGGAGTCGCTCTACG){s<=1}(?P<row_idx>.{9})(?P<DBU1>CCATACGAGCACATTACGGG){s<=2}(?P<bc1>.+)(?P<DBU2>CTAACTCGCATACCTCTGATAAC){s<=2}(?P<lox2272>.+)(?P<DBlox2272>TGTCAGCACTCTGTCAAAATAGATCG){s<=2}....(?P<DBD1>TCGATAGGTGCGTGTGAAGG){s<=2}(?P<bc2>.+)(?P<DBD2>CCTCAGTCGCTCAGTCAAG){s<=2}(?P<col_idx>.{9})(?P<PS2>CAGGACCTAAAGAATCCCATCC){s<=1}(?P<plate2>.{9}).*$
READ_FLASH_STRUCTURE_2=^.+(?P<plate1>.{9})(?P<PS1>TAACTTACGGAGTCGCTCTACG){s<=1}(?P<row_idx>.{9})(?P<DBloxP>TTAAGTAAGTAATGGCAAACAGCTATTAT){s<=2}(?P<loxp>.+){s<=2}(?P<DBU1>CCATACGAGCACATTACGGG){s<=2}(?P<bc1>.+)(?P<DBU2>CTAACTCGCATACCTCTGATAAC){s<=2}(?P<lox2272>.+){s<=2}(?P<DBlox2272>TGTCAGCACTCTGTCAAAATAGATCG){s<=2}(?P<col_idx>.{9})(?P<PS2>CAGGACCTAAAGAATCCCATCC){s<=1}(?P<plate2>.{9}).*$
READ_FLASH_STRUCTURE_3=^.+(?P<plate1>.{9})(?P<PS1>TAACTTACGGAGTCGCTCTACG){s<=1}(?P<row_idx>.{9})(?P<ADU1>CCCTTAGAACCGAGAGTGTG){s<=2}(?P<bc1>.+)(?P<ADU2>CACTCCGTTCGTCACTCAATAA){s<=2}...(?P<ADloxP>CAGCACTCTGTCAAAATAGATCGGAAA){s<=2}(?P<loxp>.+){s<=2}(?P<ADD1>CTCCAGGGTTAGGCAGATG){s<=2}(?P<bc2>.+)(?P<ADD2>CAATCGCACTATCCCGCTG){s<=2}(?P<col_idx>.{9})(?P<PS2>CAGGACCTAAAGAATCCCATCC){s<=1}(?P<plate2>.{9}).*$
READ_FLASH_STRUCTURE_4=^.+(?P<plate1>.{9})(?P<PS1>TAACTTACGGAGTCGCTCTACG){s<=1}(?P<row_idx>.{9})(?P<ADloxP>CAGCACTCTGTCAAAATAGATCGG){s<=2}...(?P<loxp>.+){s<=2}(?P<ADD1>CTCCAGGGTTAGGCAGATG){s<=2}(?P<bc2>.+)(?P<ADD2>CAATCGCACTATCCCGCTG){s<=2}(?P<lox2272>.+){s<=2}(?P<ADlox2272>TTATGGGTATTATGGGTTAGGGATAACAGGGTAATAAGCTCTAAGTAAGTAACGGCCGCC){s<=3}(?P<col_idx>.{9})(?P<PS2>CAGGACCTAAAGAATCCCATCC){s<=1}(?P<plate2>.{9}).*$

# Sequence filtering and correction for each segment
plate1.filtered.corrected=QUALITY_FILTER(source:plate1,min_nucleotide_Q-score:10,min_avg_Q-score:20) >> M2A_CORRECTION(levenshtein_distance:1, path:/path/to/example-dataset/allowlists/rcppcr_illumina_idx.txt)
plate2.filtered.corrected =QUALITY_FILTER(source:plate2,min_nucleotide_Q-score:10,min_avg_Q-score:20) >> M2A_CORRECTION(levenshtein_distance:1, path:/path/to/example-dataset/allowlists/rcppcr_illumina_idx.txt)
row_idx.filtered.corrected =QUALITY_FILTER(source:row_idx,min_nucleotide_Q-score:10,min_avg_Q-score:20) >> M2A_CORRECTION(levenshtein_distance:1, path:/path/to/example-dataset/allowlists/rcppcr_illumina_idx.txt)
col_idx.filtered.corrected =QUALITY_FILTER(source:col_idx,min_nucleotide_Q-score:10,min_avg_Q-score:20) >> M2A_CORRECTION(levenshtein_distance:1, path:/path/to/example-dataset/allowlists/rcppcr_illumina_idx.txt)

bc1.filtered=QUALITY_FILTER(source:bc1,min_nucleotide_Q-score:10,min_avg_Q-score:20)
bc2.filtered=QUALITY_FILTER(source:bc2,min_nucleotide_Q-score:10,min_avg_Q-score:20)
lox2272.filtered=QUALITY_FILTER(source:lox2272,min_nucleotide_Q-score:10,min_avg_Q-score:20)
DBlox2272.filtered=QUALITY_FILTER(source:DBlox2272,min_nucleotide_Q-score:10,min_avg_Q-score:20)
DBloxP.filtered=QUALITY_FILTER(source:DBloxP,min_nucleotide_Q-score:10,min_avg_Q-score:20)
ADloxP.filtered=QUALITY_FILTER(source:ADloxP,min_nucleotide_Q-score:10,min_avg_Q-score:20)
loxp.filtered=QUALITY_FILTER(source:loxp,min_nucleotide_Q-score:10,min_avg_Q-score:20)
ADlox2272.filtered=QUALITY_FILTER(source:ADlox2272,min_nucleotide_Q-score:10,min_avg_Q-score:20)

[demultiplex]
# Clean segment names separating the file.
KEY=plate1.filtered.corrected,plate2.filtered.corrected,row_idx.filtered.corrected,col_idx.filtered.corrected

# Export data format, FASTQ or TSV
FORMAT=FASTQ

# Output read structure if FORMAT=fastq
READ1_STRUCTURE=bc1.filtered+bc2.filtered
READ2_STRUCTURE=lox2272.filtered+loxp.filtered,ADlox2272.filt+"AAAAAAAAA"
INDEX1_STRUCTURE=
INDEX2_STRUCTURE=

# Header annotation
READ1_HEADER_ADDITION=DBlox2272.filtered,DBloxP.filtered
READ2_HEADER_ADDITION=ADloxP.filtered
```

Execution

```
Interstellar -conf /path/to/ex6.conf
```

