import os
import sys
from os.path import join
import pandas as pd

from utils.message import message_success,message_error

# relative path
snake_dir = workflow.basedir # all configfile, scripts, restructuretext, ens are relative to snakefile (this file)
working_dir = os.getcwd() # input and output path are relative to current working directory


# sub directory
config["dbDir"] = join(working_dir, config["results_name"], 'ngsdb_data')
config["genomeFasta"] = config["genomeFasta_path"]
config["genomeAnno"] = config["genomeAnno_path"]
config['exp_data'] = config['exp_path']

ngsdb_code_dir = join(os.path.dirname(snake_dir), 'ngsdb')
django_dir = join(working_dir, config['results_name'], 'ngsdb_code')

# 0. cp ngsdb to working directory

if not os.path.exists(django_dir):
    os.makedirs(django_dir, exist_ok=True)
    os.system('cp -r {}/* {}'.format(ngsdb_code_dir, django_dir))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------- #
# detail parameters in pipe #

# blast
# 1. expression matrix database create
exp_db_outdir = join(config["dbDir"], "exp")
anno_db_outdir = join(config["dbDir"], "anno")

# sqlite3
# 2. blastdb
blastdb_outdir = join(config["dbDir"], "blastdb")

# 3. gffutils
gffdb_outdir = join(config["dbDir"], "gff_sqlite3")

# 4. genomebrowse
gbrowse_outdir = join(config["dbDir"], "gbrowse")
annotation_gbrowse_outdir = join(gbrowse_outdir, 'annotation')

# 5. migration
migration_outdir = join(config["dbDir"], "migration")

# 6. addscript
addscript_outdir = join(config["dbDir"], "addscript")

rule all:
    input:
        # 1. exp database create #
        sqlite3_exp                     = join(exp_db_outdir, "exp.sqlite3"),
        exp_django_model                = join(django_dir, "geneExpAtlas", "models.py"),

        # 2. makeblastdb
        merged_db                       = join(blastdb_outdir, "merged_blastdb.ok"),

        # 3. gff database create
        gffdb                           = join(gffdb_outdir, "gtf.sqlite3"),
        gffdjango_model                 = join(django_dir, "geneAnno", "models.py"),

        # 4. genome browse
        gtfgzip                         = join(annotation_gbrowse_outdir, 'annotation.sorted.bgzip'),
        genome                          = join(annotation_gbrowse_outdir, "genome.fa"),

        # 5. migration
        makemigration                   = join(migration_outdir, "migration.ok"),

        # 6. add script to wooey
        addscript                       = join(addscript_outdir, "addscript.ok"),
        

onsuccess:
    print(message_success)
    shell("python {}/scripts/sendmail.py -r {} -t {} -l {}".format(snake_dir, config['email_addr'], "success", "{log}"))

onerror:
    print(message_error)
    shell("python {}/scripts/sendmail.py -r {} -t {} -l {}".format(snake_dir, config['email_addr'], "error", "{log}"))

include: join("rules", "8.db_generate_of_exp.Snakefile.py")
include: join("rules", "8.db_generate_of_gff.Snakefile.py")
include: join("rules", "8.db_generate_of_blastdb.Snakefile.py")
include: join("rules", "8.db_generate_of_genomebrowser.Snakefile.py")
include: join("rules", "8.db_migration.Snakefile.py")
include: join("rules", "8.db_generate_addscript.Snakefile.py")
