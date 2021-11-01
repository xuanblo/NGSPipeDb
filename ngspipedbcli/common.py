#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
define all constant variables
'''

import os
import sys
import subprocess
import time
import click
from rich import print
# colorama
from rich.console import Console
from rich.progress import track
from rich.table import Column,Table
from rich.columns import Columns
from rich.markdown import Markdown
import pandas as pd

console = Console()

#__all__ = ['modules_dict']

__version__ = "0.0.8"
ngspipedb_dir = os.path.dirname(__file__)
ngspipe_dir = os.path.join(ngspipedb_dir, 'ngspipe')
ngsdb_dir = os.path.join(ngspipedb_dir, 'ngsdb')

#dataurl = 'http://www.liu-lab.com'
dataurl = 'http://crisprlnc.org/static/download' # /home/crisprlnc/crisprlnc_website/static/download

pipes_dict = {
    'ngspipe-rnaseq-basic': {
        'snakefile': os.path.join(ngspipe_dir, '1.1.rnaseq_analysis_reference_basic.Snakefile.py'),
        'configfile': os.path.join(ngspipe_dir, 'config/rnaseq.config.yaml'),
        'env_path': os.path.join(ngspipe_dir, 'envs/requirements_rnaseq.yaml'),
        'env_name': 'ngspipe-rnaseq-basic',
        'testdata': 'testdata-ngspipe-rnaseq-basic.tar.gz',
        'database': ['eggnog.tar.gz', 'gene_ontology.tar.gz'],
        'packenv': ['ngspipe-rnaseq-basic_linux.tar.gz', 'ngspipe-rnaseq-basic_osx.tar.gz'],
        'steps': {
            'sampling_reads': [0, 'sample_path', 'rawreads_dir'],
            'rawreads_qc': [1],
            'mapping': [2, 'genomeAnno_path', 'genomeFasta_path'],
            'transcript_assembly': [3],
            'quantification': [4],
            'differential_expression': [5, 'condition_path'],
            'protein_annotation': [6, 'database_eggnog_dir', 'database_gene_ontology_path'],
            'enrich': [7],
            'network': [8],
            'all': [9]
        },
    },
    'ngspipe-tnt': {
        'snakefile': os.path.join(ngspipe_dir, 'medicago_tnt1_itis.smk'),
        'configfile': os.path.join(ngspipe_dir, 'config/tnt.config.yaml'),
        'env_path': os.path.join(ngspipe_dir, 'envs/requirements_tnt.yaml'),
        'env_name': 'ngspipe-tnt',
        'testdata': 'testdata-ngspipe-rnaseq-basic.tar.gz',
        'database': ['', ''],
        'packenv': ['', ''],
        'steps': {
            'sampling_reads': [0, 'sample_path', 'rawreads_dir'],
            'rawreads_qc': [1],
            'tnt_merge': [2, 'genomeAnno_path', 'genomeFasta_path'],
            'all': [3]
        },
    }, 
    'ngspipe-rnaseq-lncRNA': {
        'snakefile': os.path.join(ngspipe_dir, '1.3.rnaseq_analysis_reference_novelrna.Snakefile.py'),
        'configfile': os.path.join(ngspipe_dir, 'config/lncRNA.config.yaml'),
        'env_path': os.path.join(ngspipe_dir, 'envs/requirements_rnaseq.yaml'),
        'env_name': 'ngspipe-rnaseq-basic',
        'testdata': 'testdata-ngspipe-rnaseq-basic.tar.gz',
        'database': ['eggnog.tar.gz', 'kegg.tar.gz'],
        'packenv': ['ngspipe-rnaseq-basic_linux.tar.gz', 'ngspipe-rnaseq-basic_osx.tar.gz'],
    }, 
    'ngspipe-rnaseq-trinity': {
        'snakefile': os.path.join(ngspipe_dir, '1.4.rnaseq_analysis_denovo_trinity.Snakefile.py'),
        'configfile': os.path.join(ngspipe_dir, 'config/trinity.config.yaml'),
        'env_path': os.path.join(ngspipe_dir, 'envs/requirements_ngspipe_trinity.yaml'),
        'env_name': 'ngspipe-rnaseq-trinity',
        'testdata': 'testdata-ngspipe-rnaseq-basic.tar.gz',
        'database': ['eggnog.tar.gz', 'kegg.tar.gz'],
        'packenv': ['ngspipe-rnaseq-trinity_linux.tar.gz', 'ngspipe-rnaseq-trinity_osx.tar.gz'],
    }, 
    'ngspipe-chipseq': {
        'snakefile': os.path.join(ngspipe_dir, '3.1.chipseq.Snakefile.py'),
        'configfile': os.path.join(ngspipe_dir, 'config/chipseq.config.yaml'),
        'env_path': os.path.join(ngspipe_dir, 'envs/requirements_chipseq.yaml'),
        'env_name': 'ngspipe-chipseq',
        'testdata': 'testdata-ngspipe-rnaseq-basic.tar.gz',
        'database': [],
        'packenv': ['ngspipe-chipseq_linux.tar.gz', 'ngspipe-chipseq_osx.tar.gz'],
    }, 
    'ngspipe-resequencing': {
        'snakefile': os.path.join(ngspipe_dir, '2.1.resequencing_analysis.Snakefile.py'),
        'configfile': os.path.join(ngspipe_dir, 'config/resequencing.config.yaml'),
        'env_path': os.path.join(ngspipe_dir, 'envs/requirements_resequencing.yaml'),
        'env_name': 'ngspipe-resequencing',
        'testdata': 'testdata-ngspipe-rnaseq-basic.tar.gz',
        'database': [],
        'packenv': ['ngspipe-resequencing_linux.tar.gz', 'ngspipe-resequencing_linux.osx.gz'],
    },
    'ngsdb': {
        'snakefile': os.path.join(ngspipe_dir, 'db_generate.Snakefile.py'),
        'configfile': os.path.join(ngspipe_dir, 'config/ngsdb.config.yaml'),
        'env_path': os.path.join(ngspipe_dir, 'envs/requirements_ngsdb.yaml'),
        'env_name': 'ngsdb',
        'testdata': 'preparing',
        'database': [],
        'packenv': ['ngsdb_linux.tar.gz', 'ngsdb_osx.tar.gz'],
    },
}

ngspipedb_envs_set = set([pipes_dict[i]['env_name'] for i in pipes_dict.keys()])

def print_rich_markdown_data_table():
    table = Table(show_header=True, header_style='bold magenta')
    table.add_column('Pipeline', style='dim')
    table.add_column('Testdata')
    table.add_column('Database')
    table.add_column('Environment')
    for pipe in pipes_dict.keys():
        table.add_row(
            pipe,
            #Columns([pipes_dict[pipe]['testdata']]),
            pipes_dict[pipe]['testdata'],
            #Markdown('\n'.join(['- ' + i for i in pipes_dict[pipe]['database']])),
            '\n'.join(['' + i for i in pipes_dict[pipe]['database']]),
            '\n'.join(['' + i for i in pipes_dict[pipe]['packenv']]),
        )
    console.print(table)

def current_date(format=1):
    from datetime import date
    # datetime object containing current date and time
    today = date.today()
    
    if format==1:
        # dd/mm/YY
        d1 = today.strftime("%d/%m/%Y")
        #print("d1 =", d1)
        return d1
    if format==2:
        # Textual month, day and year	
        d2 = today.strftime("%B %d, %Y")
        #print("d2 =", d2)
        return d2
    if format==3:
        # mm/dd/y
        d3 = today.strftime("%m/%d/%y")
        #print("d3 =", d3)
        return d3
    if format==4:
        # Month abbreviation, day and year	
        d4 = today.strftime("%b-%d-%Y")
        #print("d4 =", d4)
        return d4

def current_datetime():
    from datetime import datetime
    # datetime object containing current date and time
    now = datetime.now()
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    #print("date and time =", dt_string)
    return dt_string


def ngspipedb_print_command(task, command):
    '''
    input a shell commond string, output a ngspipedb style [COMMAND]
    magenta: 
    cyan: 
    yellow:
    green:
    red:
    white:
    blue:
    black:
    '''
    shell_command = '[NGSPipeDb | {datetime} | {task}] {command} '.format(datetime=current_datetime(), task=task, command=command)
    #console.log(shell_command)
    #print(shell_command)
    click.secho(shell_command, bg='', fg='magenta', err=False, bold=False)

def ngspipedb_print_rich_stderr(text):
    print(text)

def ngspipedb_print_rich_stdout(text):
    print(text)

def ngspipedb_print_default_stdout(text):
    sys.stdout.write(text+'\n')

def config_validate(config):
    '''
    check if file exists
    '''
    pass

def call_snakemake_module(module_name):
    '''
    give a module, this function will retrurn snakefile and config.yaml
    snakefile is the pipeline
    '''
    return pipes_dict[module_name]

def _cmd():
    cmd = (
        "snakemake --snakefile {snakefile} --directory {working_dir} "
        "--jobs {jobs} --rerun-incomplete "
        "--configfile '{config_file}' --nolock "
        " {profile} --use-conda {conda_prefix} {dryrun} "
        " {target_rule} "
        " {args} "
    ).format(
        snakefile=get_snakefile(),
        working_dir=working_dir,
        jobs=jobs,
        config_file=config_file,
        profile="" if (profile is None) else "--profile {}".format(profile),
        dryrun="--dryrun" if dryrun else "",
        args=" ".join(snakemake_args),
        target_rule=workflow if workflow!="None" else "",
        conda_prefix= "--conda-prefix "+os.path.join(db_dir,'conda_envs')
    )
    logging.info("Executing: %s" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logging.critical(e)
        exit(1)
