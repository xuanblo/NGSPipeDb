# NGSPipeDb: NGS Pipelines & Databases

![GitHub release (latest by date)](https://img.shields.io/github/v/release/xuanblo/NGSPipeDb)
![license GPL-3.0](https://img.shields.io/github/license/xuanblo/NGSPipeDb)
![build](https://img.shields.io/travis/com/xuanblo/NGSPipeDb)

__NGSPipeDb__ is an automated pipeline for parallel processing of huge next generation sequencing (NGS) data and database generation using [snakemake workflow](https://snakemake.readthedocs.io/en/stable/index.html) which allows for ease of use, optimal speed, and a highly modular code that can be further added onto and customized by experienced users. It can be further divided into `NGSPipe` and `NGSDb` for individual usage.

## Quick start

### Installation

1. `pip3 install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple ngspipedb`
2. use conda (developing...)
3. docker (developing...)

### Usage

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

### Example:

*This example shows you how to run RNA-Seq analysis and visualize result in a database*.

1.download test data:

  ```shell
  ngspipedb download -n ngspipe-rnaseq-basic -t testdata && tar -zxvf testdata_ngspipe-rnaseq-basic.tar.gz
  ```

2.run rnaseq analysis on test data:

  ```shell
  ngspipedb runpipe ngspipe-rnaseq-basic -n ngspipe-rnaseq-basic -d test_pipeline --genomeFasta testdata_ngspipe-rnaseq-basic/genome/chr19.fa --genomeAnno testdata_ngspipe-rnaseq-basic/genome/GRCm38.83.chr19.gtf --samplefile testdata_ngspipe-rnaseq-basic/rawdata/sample.csv --conditionfile testdata_ngspipe-rnaseq-basic/rawdata/condition.csv --rawreadsdir testdata_ngspipe-rnaseq-basic/rawdata --snaketype p --report -db -ps
  ```

3.start ngsdb server:

  ```shell
  ngspipedb rundb server -url 0.0.0.0:8000
  ```

## Tutorial

A more detailed tutorial of how to use this toolkit can be found here:: [https://xuanblo.github.io/NGSPipeDb/](https://xuanblo.github.io/NGSPipeDb/)

Some demos:
- Demo RNA-Seq results: [http://www.liu-lab.com](http://www.liu-lab.com)
- Demo report: [http://www.liu-lab.com](http://www.liu-lab.com)
- Demo database: [http://www.liu-lab.com](http://www.liu-lab.com)

## Change logs

- *[2021-9-3]* update to v0.0.13a
- *[2021-9-3]* update to v0.0.13