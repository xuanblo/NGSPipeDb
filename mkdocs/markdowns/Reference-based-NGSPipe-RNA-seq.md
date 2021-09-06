---
title: 
date: 2021-03-07 09:51:19
tags:
  - Markdown
  - rnaseq
categories: module

---

# Reference-base RNA-seq analysis use NGSPipe

#### 

!!! Info inline end
    If this is your first time using NGSPipe, then we strongly recommend that you start by running test data. If you already have experience with NGSPipe, we suggest you can go straight to the custom data section.

Reference genome-based - an assembled genome exists for a species for which an RNAseq experiment is performed. It allows reads to be aligned against the reference genome and significantly improves our ability to reconstruct transcripts.

---

??? note "A typical flow of transcriptome analysis with reference is shown in the figure below"
    ![img](../imgs/RNA-seq-analysis-flow-chart-An-example-RNA-seq-analysis-workflow-is-depicted-for-a_W640.jpeg)


## Quick Start - RNA-Seq analysis on test data <a name="QuickStarted"></a>

### 1. Download test files <a name="Testdata"></a>

NGSPipe is dependent on reference files and raw sequence reads which can be found in [http://www.liu-lab.com/ngspipedb/testdata](http://www.liu-lab.com/ngspipedb/testdata).

To download the mouse RNA-seq test data into `./test_pipeline`.

```shell
ngspipedb download -n ngspipe-rnaseq-basic -t testdata -d test_pipeline
cd test_pipeline
tar -zxvf testdata-ngspipe-rnaseq-basic.tar.gz
cd ..
# -n pipeline name
# -t data type
# ngspipedb download -h for help
```

Make sure you have the following directory structure by command `tree test_pipeline`:

    testdata_ngspipe-rnaseq-basic
    ├── genome
    │   ├── GRCm38.83.chr19.gtf
    │   └── chr19.fa
    └── rawdata
        ├── condition.csv
        ├── control-0_R1.fq.gz
        ├── control-0_R2.fq.gz
        ├── control-1_R1.fq.gz
        ├── control-1_R2.fq.gz
        ├── control-2_R1.fq.gz
        ├── control-2_R2.fq.gz
        ├── sample.csv
        ├── treated-0_R1.fq.gz
        ├── treated-0_R2.fq.gz
        ├── treated-1_R1.fq.gz
        ├── treated-1_R2.fq.gz
        ├── treated-2_R1.fq.gz
        └── treated-2_R2.fq.gz

    2 directories, 16 files
!!! warning
    The test data is only used to verify that the analytical process is working properly and the analysis results do not have a biological significance.

### 2. Run RNA-seq analysis on test data <a name="RunTest"></a>

We provied a basic reference-based RNA-seq workflow for users to take a glance of NGSPipe. This workflow contains 7 steps:

    1. sampling data (choose part of your data)
    2. raw reads qc
    3. junction align to genome
    4. transcript assembly
    5. gene quantification
    6. statistic
    7. differential gene analysis

You can do RNA-seq analysis by just one simply command.

```shell
ngspipedb runpipe ngspipe-rnaseq-basic -n ngspipe-rnaseq-basic -d test_pipeline --genomeFasta testdata_ngspipe-rnaseq-basic/genome/chr19.fa --genomeAnno testdata_ngspipe-rnaseq-basic/genome/GRCm38.83.chr19.gtf --samplefile testdata_ngspipe-rnaseq-basic/rawdata/sample.csv --conditionfile testdata_ngspipe-rnaseq-basic/rawdata/condition.csv --rawreadsdir testdata_ngspipe-rnaseq-basic/rawdata --snaketype p --report -db -ps
```

The final data files are put in the folder `test_pipeline/ngspipe-rnaseq-basic`. Please check you result file `tree -d -L 2 test_pipeline/ngspipe-rnaseq-basic`, it may like this:

    test_pipeline/ngspipe-rnaseq-basic
    ├── database
    ├── genome
    ├── rawdata
    └── result_Sep-06-2021
        ├── ngsdb_code
        │   ├── __pycache__
        │   ├── blastplus
        │   ├── geneAnno
        │   ├── geneDetail
        │   ├── geneExpAtlas
        │   ├── home
        │   ├── igv
        │   ├── media
        │   ├── ngsdb
        │   ├── search
        │   ├── tools
        │   └── wooey
        ├── ngsdb_data
        │   ├── addscript
        │   ├── blastdb
        │   ├── exp
        │   ├── gbrowse
        │   ├── gff_sqlite3
        │   └── migration
        ├── ngspipe_result
        │   ├── diff
        │   ├── mapping
        │   ├── quantify
        │   ├── rawReads_qc
        │   ├── sampling_data
        │   └── statistic
        └── report
            ├── 1.pipeline
            ├── 2.rawreads_stat
            ├── 3.cleanreads_stat
            ├── 4.mapping_stat
            └── 5.exp_stat

    37 directories

!!! Note
    If you encounter any problem in this step, please turn to `TroubleShooting` for help.

---