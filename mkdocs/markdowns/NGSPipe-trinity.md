# denovo RNA-Seq analysis

## step-by-step run RNA-Seq analysis on testdata

### 1. pre-prepare download source code

1. [Install wget and git](../NGSPipe-RNA-seq/#BasicLinux)
2. [Install Miniconda3](../NGSPipe-RNA-seq/#Miniconda)
3. [Download NGSPipeDb source code](../NGSPipe-RNA-seq/#NGSPipeDbSource)

Modify the project name and enter the project directory.

```shell
mv NGSPipeDb species_sample_transcript_analysis_by_NGSPipeDb
cd species_sample_transcript_analysis_by_NGSPipeDb
```

### 2. create conda envirenment

conda bioconda lastest trinity version on macos is vdate.2011_11_26; on linux is v2.12.0

建议还是在linux上运行trinity

```shell
mamba create -n ngspipe-trinity python=3.8 -c conda-forge -y
mamba env update -n ngspipe-trinity --file ngspipe/envs/requirements_ngspipe_trinity.yaml --prune
conda activate ngspipe-trinity
```

### 3. download testdata

```
mkdir -p testdata && cd testdata
wget http://www.liu-lab.com/ngspipedb/testdata/control_R1.fq.gz
wget http://www.liu-lab.com/ngspipedb/testdata/treated_R1.fq.gz
wget http://www.liu-lab.com/ngspipedb/testdata/control_R2.fq.gz
wget http://www.liu-lab.com/ngspipedb/testdata/treated_R2.fq.gz
wget http://www.liu-lab.com/ngspipedb/testdata/chr19.fa.gz
gunzip chr19.fa.gz
wget http://www.liu-lab.com/ngspipedb/testdata/samples_resequencing.xls -O samples_chipseq.xls
cd ..
```

### 4. run snakemake

```
snakemake -s ngspipe/1.4.rnaseq_analysis_denovo_trinity.Snakefile.py --configfile ngspipe/config/trinity.config.yaml -p -j 1 -n
```
