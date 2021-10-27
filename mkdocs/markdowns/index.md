# NGSPipeDb: NGS Pipelines & Databases

__NGSPipeDb__ is an automated pipeline for parallel processing of huge next generation sequencing (NGS) data and database generation using [snakemake workflow](https://snakemake.readthedocs.io/en/stable/index.html) which allows for ease of use, optimal speed, and a highly modular code that can be further added onto and customized by experienced users.

<figure markdown> 
  ![Dummy image](imgs/fig1.png){ width="800" }
  <figcaption>Overview of NGSPipeDb</figcaption>
</figure>

## Quick start

### 0. required

Although included in this section are step-by-step instructions, it is assumed that the user has a basic understanding of the [nix command line interface](https://en.wikipedia.org/wiki/Command-line_interface). Also, it would be better if the user has basic knowledge about [snakemake](https://snakemake.readthedocs.io/en/stable/), [conda](https://docs.conda.io/en/latest/) and [best practice RNA sequence analysis](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8), but it is not required. You can also find some easy-to-learn matierals in our "Learning materials" page, for example [linux & shell](../linux) and [RNASeq background](../ngs#rnaseq) for beginers.

**Install wget and git** <a name="BasicLinux"></a>

To get some of the required software packages, we will use the command line tools called [wget](http://www.gnu.org/software/wget/) and [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).  *wget* is a popular tool for downloading things off of the internet.  *git* is a distributed version control system which we will use to checkout the NGSPipeDb code.

!!! note
    These tools are already pre-installed in most systems, but if you are unsure whether or not you have *wget* enter `wget` and if the return is `wget: command not found`, then you will have to install *wget*.  Do likewise for *git*.

**Install Miniconda3** <a name="Miniconda"></a>

NGSPipeDb relies on the conda package manager for installation and dependency resolution, so you will need to [install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) first.

We will be using the [Miniconda3](http://conda.pydata.org/miniconda.html) package management system (aka CONDA) to manage all of the software packages that NGSPipe is dependent on. 

Use following commands to retrieve and then run the Minicoda3 installation script:

1.download miniconda3

=== "Linux & WSL"

    ```shell
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    ```

=== "MacOSX"

    ```shell
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    ```

2.install miniconda3

=== "Linux & WSL"

    ```shell
    bash Miniconda3-latest-Linux-x86_64.sh
    ```

=== "MacOSX"

    ```shell
    bash Miniconda3-latest-MacOSX-x86_64.sh
    ```

!!! important
    While running the installation script, follow the commands listed on screen, and press the _enter_ key to scroll. Make sure to answer `yes` when asked if you want to prepend Miniconda3 to PATH. After that, close your terminal, open a new one and you should now have Conda working! You could, alternatively, run `#!shell source ~/.bashrc` to initiate conda.

3.Test if conda is ready to work by entering: `conda update conda`. Press `y` to confirm the conda updates.

4.Finally, `conda install mamba -c conda-forge`.

!!! Note
    Mamba is a reimplementation of the conda package manager in C++, the fast conda-alternative. Mamba is recommended but not necessary.

!!! Info
    You will only have to install Minicoda3 once. If you face any conda problem, please [learn](../conda) more about it.

### Installation

1. Install from pipi
  ```shell
  conda install pip
  pip3 install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple ngspipedb
  ```

2. Install from conda
```shell
conda install
```

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

**run resequencing analysis, generate report, and build RNA-Seq database**

## Tutorial

A more detailed tutorial of how to use this toolkit can be found here:: [https://xuanblo.github.io/NGSPipeDb/](https://xuanblo.github.io/NGSPipeDb/)

Some demos:
- Demo RNA-Seq report: [http://rnaseq-report.liu-lab.com](http://rnaseq-report.liu-lab.com)
- Demo RNA-Seq database: [http://ngsdb-rnaseq.liu-lab.com](http://ngsdb-rnaseq.liu-lab.com)

## Change logs

- https://xuanblo.github.io/NGSPipeDb/changelog/