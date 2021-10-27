# Medicago tnt1 insertion

> - ITIS: tool to Identify TE Insertion Sites in genome
> - https://github.com/Chuan-Jiang/ITIS
> - https://pubmed.ncbi.nlm.nih.gov/25887332/

Example:
```
python -m ngspipedbcli runpipe ngspipe-tnt-medicago -n ngspipe-tnt --resultdirname result -d ../ --genomeFasta ../genome/scf003.fa --genomeAnno ../genome/scf003.sorted.gff --samplefile ../rawdata/sample.csv --rawreadsdir ../rawdata --snaketype p --reads_prefix _R{}_paired.fastq.gz -j 20
```