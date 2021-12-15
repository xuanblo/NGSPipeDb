# galaxy #

# rnaseq

## test command line get parameters
python -m ngspipedbcli runpipe ngspipe-rnaseq-basic -d ../test_pipeline -n ngspipe-rnaseq-basic --resultdirname result --genomeFasta ../testdata_ngspipe-rnaseq-basic/genome/chr19.fa --genomeAnno ../testdata_ngspipe-rnaseq-basic/genome/GRCm38.83.chr19.gtf --samplefile ../testdata_ngspipe-rnaseq-basic/rawdata/sample.csv --conditionfile ../testdata_ngspipe-rnaseq-basic/rawdata/condition.csv --rawreadsdir ../testdata_ngspipe-rnaseq-basic/rawdata --target differential_expression --snaketype p -j 1
## test get parameters from configfile
python -m ngspipedbcli runpipe ngspipe-rnaseq-basic -d ../test_pipeline -n ngspipe-rnaseq-basic -c ../test_pipeline/ngspipe-rnaseq-basic/ngspipe_config.yaml --resultdirname result --snaketype p -j 1 -r

# medicago tnt

## parameters
python -m ngspipedbcli runpipe ngspipe-tnt-medicago -n ngspipe-tnt --resultdirname result -d ../test_pipeline --genomeFasta ../testdata_ngspipe-rnaseq-basic/genome/chr19.fa --genomeAnno ../testdata_ngspipe-rnaseq-basic/genome/GRCm38.83.chr19.gtf --exogenous_seq ../testdata_ngspipe-rnaseq-basic/database/exogenous/medicago_tnt1.fa --samplefile ../testdata_ngspipe-rnaseq-basic/rawdata/sample.csv --rawreadsdir ../testdata_ngspipe-rnaseq-basic/rawdata --snaketype np --readsprefix _R{}.fq.gz -j 20
## config
python -m ngspipedbcli startproject ngspipe-tnt-medicago -n ngspipe-tnt -d ../test_pipeline
vi ../test_pipelines/ngspipe-tnt-medicago/
python -m ngspipedbcli runpipe ngspipe-tnt-medicago -n ngspipe-tnt -d ../test_pipeline -c ../test_pipeline/ngspipe-tnt-medicago/ngspipe_config.yaml --snaketype np --readsprefix _R{}.fq.gz -j 1

## resequencing
python -m ngspipedbcli runpipe ngspipe-resequencing -n ngspipe-resequencing -d ../test_pipeline --genomeFasta ../testdata_ngspipe-rnaseq-basic/genome/chr19.fa --genomeAnno ../testdata_ngspipe-rnaseq-basic/genome/GRCm38.83.chr19.gtf --samplefile ../testdata_ngspipe-rnaseq-basic/rawdata/sample.csv --rawreadsdir ../testdata_ngspipe-rnaseq-basic/rawdata --snaketype np --readsprefix _R{}.fq.gz -j 1


# trinity

python -m ngspipedbcli runpipe ngspipe-rnaseq-trinity -d ../test_pipeline -n ngspipe-rnaseq-trinity --resultdirname result_Sep-18-2021 --samplefile ../testdata_ngspipe-rnaseq-basic/rawdata/sample.csv --conditionfile ../testdata_ngspipe-rnaseq-basic/rawdata/condition.csv --rawreadsdir ../testdata_ngspipe-rnaseq-basic/rawdata --snaketype np -ps

# lncRNA

python -m ngspipedbcli runpipe ngspipe-rnaseq-lncRNA -d ../test_pipeline -n ngspipe-rnaseq-lncRNA --resultdirname result_Sep-18-2021 --genomeFasta ../testdata_ngspipe-rnaseq-basic/genome/chr19.fa --genomeAnno ../testdata_ngspipe-rnaseq-basic/genome/GRCm38.83.chr19.gtf --samplefile ../testdata_ngspipe-rnaseq-basic/rawdata/sample.csv --conditionfile ../testdata_ngspipe-rnaseq-basic/rawdata/condition.csv --rawreadsdir ../testdata_ngspipe-rnaseq-basic/rawdata --snaketype np -ps

# ngsdb

python -m ngspipedbcli rundb build ngsdb --genomeFasta ../testdata_ngspipe-rnaseq-basic/genome/chr19.fa --genomeAnno ../testdata_ngspipe-rnaseq-basic/genome/GRCm38.83.chr19.gtf -exp ../testdata_ngspipe-rnaseq-basic/db/gene_fpkm_all_samples.tsv --resultdirname result -d ../test_pipeline -ps
python -m ngspipedbcli rundb build ngsdb --genomeFasta ../testdata_ngspipe-rnaseq-basic/genome/chr19.fa --genomeAnno ../testdata_ngspipe-rnaseq-basic/genome/GRCm38.83.chr19.gtf -exp ../testdata_ngspipe-rnaseq-basic/db/gene_fpkm_all_samples.tsv --resultdirname result -d ../test_pipeline --snaketype p
python -m ngspipedbcli rundb serve -m ../test_pipeline/ngsdb/result/ngsdb_code/manage.py -up 127.0.0.1:8000 -ps

# env
python -m ngspipedbcli env update --pipename ngspipe-rnaseq-basic -n ngspipe-rnaseq-basic

# wsgi

https://docs.djangoproject.com/zh-hans/2.2/howto/deployment/wsgi/uwsgi/
conda install uwsgi

uwsgi --reload /tmp/ngsdb_demo.pid

mkdir commonstatic
(base) [zhangxuan@chengdu ngsdb_code]$ python manage.py collectstatic

datatable 报ajax错误，在ini中增加buffer

https://www.cnblogs.com/ryxiong-blog/p/12894077.html

工具在nginx和uwsgi上运行不了