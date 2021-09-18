# NGSPipeDb: NGS Pipelines & Databases

__NGSPipeDb__ is an automated pipeline for parallel processing of huge next generation sequencing (NGS) data and database generation using [snakemake workflow](https://snakemake.readthedocs.io/en/stable/index.html) which allows for ease of use, optimal speed, and a highly modular code that can be further added onto and customized by experienced users.

## Quick start

### Installation

1. Install from pipi: `pip3 install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple ngspipedb`
2. Install from conda: `use conda (developing...)`

### NGSPipeDb command line interface

`ngspipedb -h`

```yaml
Usage: ngspipedb [OPTIONS] COMMAND [ARGS]...

  ngspipedb is a snakemake-based tool for reproducible next generation
  sequencing (NGS) data analysis and interactive web application auto-build.

  Example:
  ngspipedb env create -n ngspipe-rnaseq-basic
  ngspipedb download -n ngspipe-rnaseq-basic -t testdata
  ngspipedb startproject myprojectname -n ngspipe-rnaseq-basic
  ngspipedb runpipe myprojectname ngspipe-rnaseq-basic --report -db
  ngspipedb rundb serve test_pipeline/ngspipe-rnaseq-basic/result_Sep-06-2021/ngsdb_code/manage.py -up 0.0.0.0:8000

  A more detailed tutorial of how to use this toolkit can be found here:
  https://xuanblo.github.io/NGSPipeDb/

Options:
  --version   Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  download      data retrive related commands
  env           ngspipedb environment related commands
  rundb         generate a database related commands
  runpipe       run ngspipe
  startproject  Creates a ngspipedb project directory structure for the given
                project name in the current directory or optionally in the
                given directory.
```

### Usage example: run RNA-Seq analysis, generate report, and build RNA-Seq database

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
  ngspipedb rundb server -url 0.0.0.0:8000
  ```

## Tutorial

A more detailed tutorial of how to use this toolkit can be found here:: [https://xuanblo.github.io/NGSPipeDb/](https://xuanblo.github.io/NGSPipeDb/)

Some demos:
- Demo RNA-Seq results: [http://rnaseq-result.liu-lab.com](http://rnaseq-result.liu-lab.com)
- Demo report: [http://rnaseq-report.liu-lab.com](http://rnaseq-report.liu-lab.com)
- Demo RNA-Seq database: [http://ngsdb-rnaseq.liu-lab.com](http://ngsdb-rnaseq.liu-lab.com)

## Change logs

- *[2021-9-18]* update to v0.0.16
- *[2021-9-3]* update to v0.0.13a
- *[2021-9-3]* update to v0.0.13