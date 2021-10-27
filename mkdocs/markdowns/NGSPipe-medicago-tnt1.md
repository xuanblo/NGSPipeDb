# Medicago tnt1 insertion

> - ITIS: tool to Identify TE Insertion Sites in genome
> - https://github.com/Chuan-Jiang/ITIS
> - https://pubmed.ncbi.nlm.nih.gov/25887332/

Example:
```
python -m ngspipedbcli runpipe ngspipe-tnt-medicago -n ngspipe-tnt --resultdirname result -d ../ --genomeFasta ../genome/scf003.fa --genomeAnno ../genome/scf003.sorted.gff --samplefile ../rawdata/sample.csv --rawreadsdir ../rawdata --snaketype p --reads_prefix _R{}_paired.fastq.gz -j 20
```

- `ngspipe-tnt-medicago` your project name
- `-n ngspipe-tnt` pipeline name
- `--genomeFasta ../genome/scf003.fa` give a genome fasta file path, see file format [fasta](https://en.wikipedia.org/wiki/FASTA_format)
- `--genomeAnno ../genome/scf003.sorted.gff` give a genome annotaion file path [gtf](https://genome.ucsc.edu/FAQ/FAQformat.html#format4)/[gff](https://genome.ucsc.edu/FAQ/FAQformat.html#format3)
- `--samplefile ../rawdata/sample.csv` give a sample file path, which has one column

samplefile:
```yaml
XTD-6
```

configfile:
```yaml
#---------------------------
# medicago tnt1
#---------------------------

## 1.reference ##
genomeAnno_path: /home/zhangxuan/Work/Project/hehua/test_ngspipedb/genome/scf003.sorted.gff # gene annotation file, can be gtf or gff format
genomeFasta_path: /home/zhangxuan/Work/Project/hehua/test_ngspipedb/genome/scf003.fa # genome sequence, fasta format

## 2.raw reads data ##
sample_path: /home/zhangxuan/Work/Project/hehua/test_ngspipedb/rawdata/sample.csv # sample file
rawreads_dir: /home/zhangxuan/Work/Project/hehua/test_ngspipedb/rawdata # sample file directory
read1Suffix: _R1_paired.fastq.gz # fastq suffix, read1
read2Suffix: _R2_paired.fastq.gz

## 3.output directory ##
results_name: result

## 4.notice ##
# if the string is 'nobody', ngspipe will not send email
# modify 'noboby' to 'xxx@qq.com' or 'xxx@qq.com,yyy@qq.com' to send email
email_addr: nobody

# choose where to stop your pipeline
target: all #

#----------------------------------
# Configuration for sampling data
#----------------------------------

# for test the pipe, you can choose to the part of the input file
# which sampling method do you want to use?
sampling_method: links # links or head or tail or seqkit_number or seqkit_proportion
# Default is links (ues the whole data of sample); head (use first sampling_range line in every sample),tail (use last sampling_range line in every sample); seqkit_number (number of reads); seqkit_proportion (percentage of reads)
# and how many reads file line or reads number or reads proportion do you want to use?
sampling_value: 80000 # for head and tail, this value is line number; for number, this value is reads number; for proportion, this value is percentage

samples_num: all # all or interger
# Default is all (use all samples), give a sample number, must less than real sample number, for example 6

#----------------------------------
# Configuration for Quality Control
#----------------------------------
# which qc method do you want to use?
qc_method: trim-galore # trim-galore or trimmomatic or fastqc
```