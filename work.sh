# galaxy #

# rnaseq

python -m ngspipedbcli runpipe ngspipe-rnaseq-basic -d ../test_pipeline -n ngspipe-rnaseq-basic --resultdirname result_Sep-06-2021 --genomeFasta ../testdata_ngspipe-rnaseq-basic/genome/chr19.fa --genomeAnno ../testdata_ngspipe-rnaseq-basic/genome/GRCm38.83.chr19.gtf --samplefile ../testdata_ngspipe-rnaseq-basic/rawdata/sample.csv --conditionfile ../testdata_ngspipe-rnaseq-basic/rawdata/condition.csv --rawreadsdir ../testdata_ngspipe-rnaseq-basic/rawdata --snaketype np


# trinity

python -m ngspipedbcli runpipe ngspipe-rnaseq-trinity -d ../test_pipeline -n ngspipe-rnaseq-trinity --resultdirname result_Sep-18-2021 --samplefile ../testdata_ngspipe-rnaseq-basic/rawdata/sample.csv --conditionfile ../testdata_ngspipe-rnaseq-basic/rawdata/condition.csv --rawreadsdir ../testdata_ngspipe-rnaseq-basic/rawdata --snaketype np -ps

# lncRNA

python -m ngspipedbcli runpipe ngspipe-rnaseq-lncRNA -d ../test_pipeline -n ngspipe-rnaseq-lncRNA --resultdirname result_Sep-18-2021 --genomeFasta ../testdata_ngspipe-rnaseq-basic/genome/chr19.fa --genomeAnno ../testdata_ngspipe-rnaseq-basic/genome/GRCm38.83.chr19.gtf --samplefile ../testdata_ngspipe-rnaseq-basic/rawdata/sample.csv --conditionfile ../testdata_ngspipe-rnaseq-basic/rawdata/condition.csv --rawreadsdir ../testdata_ngspipe-rnaseq-basic/rawdata --snaketype np -ps

# ngsdb

python -m ngspipedbcli rundb build ngsdb --genomeFasta ../testdata_ngspipe-rnaseq-basic/genome/chr19.fa --genomeAnno ../testdata_ngspipe-rnaseq-basic/genome/GRCm38.83.chr19.gtf -exp ../testdata_ngspipe-rnaseq-basic/db/gene_fpkm_all_samples.tsv --resultdirname result -d ../test_pipeline -ps
python -m ngspipedbcli rundb build ngsdb --genomeFasta ../testdata_ngspipe-rnaseq-basic/genome/chr19.fa --genomeAnno ../testdata_ngspipe-rnaseq-basic/genome/GRCm38.83.chr19.gtf -exp ../testdata_ngspipe-rnaseq-basic/db/gene_fpkm_all_samples.tsv --resultdirname result -d ../test_pipeline --snaketype p
python -m ngspipedbcli rundb serve -m ../test_pipeline/ngsdb/result/ngsdb_code/manage.py -up 127.0.0.1:8000 -ps