# Datasets

## RNA-Seq analysis

### Database

**EggNOG**  
```shell
conda activate ngspipe-rnaseq-basic
download_eggnog_data.py -y --data_dir ./
tar -zxvf eggnog.db.gz
```

**gene_ontology**
```shell
ngspipedb download -n ngspipe-rnaseq-basic -t database
tar -zxvf gene_ontology.tar.gz
```

`tree database` looks like:

    database/
    ├── eggnog
    │   ├── eggnog.db
    │   ├── eggnog_proteins.dmnd
    │   ├── eggnog.taxa.db
    │   ├── eggnog.taxa.db.traverse.pkl
    │   └── wget-log
    └── gene_ontology
        ├── go-basic.obo
        ├── go.obo
        ├── goslim_generic.obo
        └── goslim_plant.obo

    2 directories, 9 files