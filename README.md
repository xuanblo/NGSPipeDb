# NGSPipeDb: NGS Pipelines & Databases

__NGSPipeDb__ is an automated pipeline for parallel processing of huge next generation sequencing (NGS) data and database generation using [snakemake workflow](https://snakemake.readthedocs.io/en/stable/index.html) which allows for ease of use, optimal speed, and a highly modular code that can be further added onto and customized by experienced users.

## Quick start

### Required

- conda
- pip

### Installation

1. Install from pipi
  ```shell
  pip3 install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple ngspipedb
  ```

### Usage example

**run RNA-Seq analysis, generate report, and build RNA-Seq database**

step1. download test data:

  ```shell
  ngspipedb download -n ngspipe-rnaseq-basic -t testdata && tar -zxvf testdata-ngspipe-rnaseq-basic.tar.gz
  ```

step2. run rnaseq analysis on test data:

  ```shell
  ngspipedb runpipe mouse_rnaseq_analysis -n ngspipe-rnaseq-basic --genomeFasta testdata-ngspipe-rnaseq-basic/genome/chr19.fa --genomeAnno testdata-ngspipe-rnaseq-basic/genome/GRCm38.83.chr19.gtf --samplefile testdata-ngspipe-rnaseq-basic/rawdata/sample.csv --conditionfile testdata-ngspipe-rnaseq-basic/rawdata/condition.csv --rawreadsdir testdata-ngspipe-rnaseq-basic/rawdata -j 10 --report -db
  ```

step3. start ngsdb server:

  ```shell
  ngspipedb rundb serve -m mouse_rnaseq_analysis/result/ngsdb_code/manage.py -up 127.0.0.1:8000
  ```

## Tutorial

A more detailed tutorial of how to use this toolkit can be found here:: [https://xuanblo.github.io/NGSPipeDb/](https://xuanblo.github.io/NGSPipeDb/)

Some demos:
- Demo RNA-Seq report: [http://rnaseq-report.liu-lab.com](http://rnaseq-report.liu-lab.com)
- Demo RNA-Seq database: [http://ngsdb-rnaseq.liu-lab.com](http://ngsdb-rnaseq.liu-lab.com)

## Change logs

- https://xuanblo.github.io/NGSPipeDb/changelog/