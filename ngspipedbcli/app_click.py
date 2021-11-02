#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import click
from click.decorators import version_option

from ngspipedbcli.download import download_main
from ngspipedbcli.runpipe import run_ngsdb, run_ngspipe, run_ngspipedb_one_step, run_server
from ngspipedbcli.project import start_project_main
from ngspipedbcli.env import env_create, env_list, env_pack, env_remove, env_unpack, env_update
from ngspipedbcli.common import *

import requests
import json
try:
    from packaging.version import parse
except ImportError:
    from pip._vendor.packaging.version import parse
from versions_comparison import Comparison

def configure_env_group(args=None):
    '''
    TODO: @click.option('--log_level', type=click.Choice(['debug', 'info']), default='info', help='can be info debug error')
    TODO: add color to ngspipedb env list
    '''
    @click.group(name='env', short_help='ngspipedb environment related commands')
    @click.pass_context
    def cli_env(ctx, **kargs):
        """
        ngspipedb environment related commands

        Example:

        python -m ngspipedbcli env create -n ngspipe-rnaseq-basic -ps
        """
        click.echo('We are under env sub commandas')

    @cli_env.command(name='create')
    @click.pass_context
    @click.option('-n', '--name', help='ngspipedb env name', type=click.Choice(pipes_dict.keys()), required=True)
    #@click.option('-f', '--how', help='how to create environment', type=click.Choice(['conda', 'mamba', 'unpack']), default='conda')
    @click.option('-ps', '--printshell', is_flag=True, help="print ngspipedb shell commands")
    def create_cmd(ctx, **kargs):
        '''
        Create an environment based on ngspipedb task

        Example:

        python -m ngspipedbcli env create -n ngspipe-rnaseq-basic -ps
        '''
        if args:
            click.echo(ctx.get_help())
        #click.echo("ctx: {}, param".format(ctx.params))
        env_create(ctx.params)

    @cli_env.command(name='list', short_help='List the Conda environments, ngspipedb\'s environment will colored as red', epilog='真的吗')
    @click.pass_context
    def list_cmd(ctx, **kargs):
        #click.echo("ctx: {}, param".format(ctx.params))
        env_list(ctx.params)

    @cli_env.command(name='pack', help='Pack an environment')
    @click.pass_context
    @click.option('-n', '--name', help='env name', required=True)
    @click.option('-o', '--directory', type=click.Path(), help="A file name or file path", default='')
    @click.option('-p', '--platform', type=click.Choice(['osx', 'linux']), help="A file name or file path", default='linux')
    @click.option('-ps', '--printshell', is_flag=True, help="print ngspipedb shell commands")
    def pack_cmd(ctx, **kargs):
        if args:
            click.echo(ctx.get_help())
        #click.echo("ctx: {}, param".format(ctx.params))
        env_pack(ctx.params)
    
    @cli_env.command(name='remove', help='Remove an environment')
    @click.pass_context
    @click.option('-n', '--name', help='env name')
    @click.option('-p', '--path', help='env path')
    @click.option('-ps', '--printshell', is_flag=True, help="print ngspipedb shell commands")
    def pack_cmd(ctx, **kargs):
        if args:
            click.echo(ctx.get_help())
        #click.echo("ctx: {}, param".format(ctx.params))
        env_remove(ctx.params)
    
    @cli_env.command(name='unpack', help='Unpack conda environment')
    @click.pass_context
    @click.option('-f', '--file', type=click.Path(exists=True), help='packed env file', required=True)
    @click.option('-d', '--condaenvdir', help="conda env directory", default='~/miniconda3/envs')
    @click.option('-ps', '--printshell', is_flag=True, help="print ngspipedb shell commands")
    def unpack_cmd(ctx, **kargs):
        if args:
            click.echo(ctx.get_help())
        #click.echo("ctx: {}, param".format(ctx.params))
        env_unpack(ctx.params)
    
    @cli_env.command(name='update', help='update conda environment')
    @click.pass_context
    @click.option('-n', '--name', help='pipeline/env name')
    @click.option('-p', '--path', help="env path")
    @click.option('-ps', '--printshell', is_flag=True, help="print ngspipedb shell commands")
    def update_cmd(ctx, **kargs):
        if args:
            click.echo(ctx.get_help())
        #click.echo("ctx: {}, param".format(ctx.params))
        env_update(ctx.params)

    return cli_env

def configure_download_group(args=None):
    '''
    provide get/download function:
    - testdata
    - database
    TODO: @click.option('--log_level', type=click.Choice(['debug', 'info']), default='info', help='can be info debug error')
    '''
    @click.command(name='download', short_help='data retrive related commands')
    @click.pass_context
    @click.option('-l', '--list', is_flag=True, help="list all available files.")
    @click.option('-a', '--all', is_flag=True, help="download all datatypes")
    @click.option('-ps', '--printshell', is_flag=True, help="print ngspipedb shell commands")
    @click.option('-o', '--directory', type=click.Path(), default='./')
    @click.option('-n', '--pipeline', help='ngspipedb env name', type=click.Choice(pipes_dict.keys()), default='ngspipe-rnaseq-basic')
    @click.option('-p', '--platform', type=click.Choice(['osx', 'linux']), help="A file name or file path", default='linux')
    @click.option('-t', '--datatype', type=click.Choice(['env', 'testdata', 'database']), help="file types", default='testdata')
    def cli_download_cmd(ctx, **kargs):
        """
        Commands related to get testdata and database

        Example:

        python -m ngspipedbcli downloaddata -l

        python -m ngspipedbcli downloaddata -n ngspipe-rnaseq-basic -t testdata -o run_test/myproject_rnaseq_basic -ps
        """
        #click.echo('We are under env sub commandas')
        #click.echo(ctx.params)
        #click.echo(args)
        
        if args:
            click.echo(ctx.get_help())
            sys.exit()
        download_main(ctx.params)

    return cli_download_cmd

def configure_ngspipe_group(args=None):
    '''
    run ngspipe:
    - rnaseq_basic
    - rnaseq_lncRNA
    - chipseq
    TODO: check ngspipedb env if exists, auto install env
    TODO: auto create configfile in given directory
    '''
    
    @click.command(name='runpipe', short_help='run ngspipe')
    @click.argument('projectname')
    @click.pass_context
    @click.option('-n', '--pipename', help='ngspipedb env name', type=click.Choice(pipes_dict.keys()), required=True)
    @click.option('-d', '--directory', help="project directory", default='', type=click.Path())
    @click.option('-j', '--jobs', type=int, help="how many cpu to use", default=1)
    @click.option('--genomeFasta', help="genome sequence file (fasta)", default='')
    @click.option('--genomeAnno', help="genome annotation file (gff/gtf)", default='')
    @click.option('--samplefile', help="samplefile", default='')
    @click.option('--conditionfile', help="conditionfile", default='')
    @click.option('--rawreadsdir', help="raw reads directory", default='')
    @click.option('--eggnogdir', help="eggnog database directory", default='')
    @click.option('--ontologyfile', help="gene ontology database file (obo)", default='')
    @click.option('-e', '--email_addr', help="result directory name (under project directory)", default='')
    @click.option('--reads_prefix', help="reads prefix (Example: _R{}.fq.gz )", default='')
    @click.option('--resultdirname', help="result directory name under proect name directory", default='{name}_{date}'.format(name='result', date=current_date(4)))
    @click.option('--snaketype', help="`p`: print snakemake shell commands. `np`: Enable the dry run.", type=click.Choice(['np', 'p']), default='p')
    @click.option('-r', '--report', help="generate html report", is_flag=True)
    @click.option('-db', '--database', help="generate database", is_flag=True)
    @click.option('--update_env', help="install and update env auto", is_flag=True)
    @click.option('--target', help="choose where to stop your pipeline")
    @click.option('-c', '--configfile', help="config file path", type=click.Path(exists=True))
    @click.option('--otherparams', help="other snakemake params", default='')
    @click.option('-ps', '--printshell', is_flag=True, help="print ngspipedb shell commands")
    def run_pipe_cmd(ctx, **kargs):
        '''
        run a lot of pipeline

        Example:

        python -m ngspipedbcli runpipe ngspipe-rnaseq-basic -n ngspipe-rnaseq-basic -d test_pipeline --genomeFasta testdata_ngspipe-rnaseq-basic/genome/chr19.fa --genomeAnno testdata_ngspipe-rnaseq-basic/genome/GRCm38.83.chr19.gtf --samplefile testdata_ngspipe-rnaseq-basic/rawdata/sample.csv --conditionfile testdata_ngspipe-rnaseq-basic/rawdata/condition.csv --rawreadsdir testdata_ngspipe-rnaseq-basic/rawdata --snaketype p --report -db -ps

        '''
        if args:
            click.echo(ctx.get_help())
            sys.exit()
        #output = str(ctx.params)
        #click.secho(output, bg='', fg='cyan')
        run_ngspipe(ctx.params)
    
    return run_pipe_cmd

def configure_ngsdb_group(args=None):
    '''
    run ngsdb:
    - exp
    - download
    '''
    @click.group(name='rundb', short_help='generate a database related commands')
    @click.pass_context
    def cli_ngsdb(ctx, **kargs):
        '''
        Document:

        http://www.liu-lab.com
        '''
        pass
    
    @cli_ngsdb.command(name='build', short_help='interactive web application auto-build')
    @click.pass_context
    @click.argument('projectname')
    @click.option('--genomeFasta', help="genome sequence file (fasta)", default='')
    @click.option('--genomeAnno', help="genome annotation file (gff/gtf)", default='')
    @click.option('-exp', '--exp_profile', help="expression matrix, csv")
    @click.option('-c', '--configfile', help="config file path", type=click.Path(exists=True))
    @click.option('--resultdirname', help="result directory name under proect name directory", default='{name}_{date}'.format(name='result', date=current_date(4)))
    @click.option('-d', '--directory', required=True, help="working directory")
    @click.option('-e', '--email_addr', help="result directory name (under project directory)", default='')
    @click.option('-j', '--jobs', type=int, help="how many cpu to use", default=1)
    @click.option('--otherparams', help="other snakemake params", default='')
    @click.option('--snaketype', help="`p`: print snakemake shell commands. `np`: Enable the dry run.", type=click.Choice(['np', 'p']), default='p')
    @click.option('-ps', '--printshell', is_flag=True, help="print ngspipedb shell commands")
    def db_build_cmd(ctx, **kargs):
        '''
        Example:

        \b
        python -m ngspipedbcli rundb build ngspipe-rnaseq-basic -d test_pipeline --resultdirname result_Sep-04-2021 -exp test_pipeline/ngspipe-rnaseq-basic/result_Sep-04-2021/tmp_result/quantify/quantify_by_stringtie/gene_fpkm_all_samples.tsv --snaketype p
        '''
        if args:
            click.echo(ctx.get_help())
            sys.exit()
        #print(ctx.params)
        run_ngsdb(ctx.params)

    @cli_ngsdb.command(name='serve', short_help='Starts a lightweight Web server for development.')
    @click.pass_context
    @click.option('-m', '--managepy', type=click.Path(exists=True), help="manage.py path", default='./manage.py')
    @click.option('-up', '--urlport', type=str, help="web url and port", default='127.0.0.1:8000')
    @click.option('-ps', '--printshell', is_flag=True, help="print ngspipedb shell commands")
    def db_serve_cmd(ctx, **kargs):
        '''
        Example:

        \b
        python -m ngspipedbcli rundb serve -m test_pipeline/ngspipe-rnaseq-basic/result_Sep-06-2021/ngsdb_code/manage.py -up 0.0.0.0:8909
        '''
        if args:
            click.echo(ctx.get_help())
            sys.exit()
        #print(ctx.params)
        run_server(ctx.params)

    return cli_ngsdb

def configure_startproject_group(args=None):
    '''
    copy from `django-admin startproject`
    '''

    @click.command(name='startproject', short_help='Creates a ngspipedb project directory structure for the given project name in the current directory or optionally in the given directory.')
    @click.argument('projectname')
    @click.option('-n', '--pipeline', help='pipelines from ngspipedb', type=click.Choice(pipes_dict.keys()), default='ngspipe-rnaseq-basic')
    @click.pass_context
    @click.option('-d', '--directory', default='', help="project directory")
    @click.option('-ps', '--printshell', is_flag=True, help="print ngspipedb shell commands")
    def startproject_cmd(ctx, **kargs):
        '''
        Creates a ngspipedb project directory structure for the given project name in the current directory or optionally in the given directory.

        Example:

        python -m ngspipedbcli startproject myproject_rnaseq_basic -n ngspipe-rnaseq-basic -ps
        '''
        if args:
            click.echo(ctx.get_help())
        #print(ctx.params)
        start_project_main(ctx.params)

    return startproject_cmd

def configure_runserver_group(args=None):
    '''
    runserver
    '''

    @click.command(name='runserver', short_help='Starts a lightweight Web server for development.')
    @click.option('-d', '--directory', required=True, help="create project directory")
    @click.option('--dryrun', is_flag=True, help="Enable the dry run.")
    @click.pass_context
    def cli_runserver(ctx, **kargs):
        '''
        Example:

        ngspipedb startproject rnaseq-basic -d xxx
        '''
        print(ctx.params)

    return cli_runserver

def configure_info_group(args=None):
    '''
    runserver
    '''

    @click.command(name='info', short_help='see what ngspipedb can do')
    @click.option('-d', '--directory', required=True, help="create project directory")
    @click.option('--dryrun', is_flag=True, help="Enable the dry run.")
    @click.pass_context
    def cli_info(ctx, **kargs):
        '''
        Example:

        ngspipedb startproject rnaseq-basic -d xxx
        '''
        print(ctx.params)

    return cli_info

# #############################################################################################
#
# sub-parsers
#
# #############################################################################################
current_version = '0.0.24'
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.group(context_settings=CONTEXT_SETTINGS, invoke_without_command=True)
@click.version_option(version=current_version, prog_name='NGSPipeDb', message='%(prog)s, version %(version)s')
@click.pass_context
def cli(ctx, args=None):
    """
    ngspipedb is a snakemake-based tool for reproducible next generation sequencing (NGS) data analysis and interactive web application auto-build.

    \b
    Example:
    ngspipedb env create -n ngspipe-rnaseq-basic
    ngspipedb download -n ngspipe-rnaseq-basic -t testdata
    ngspipedb startproject myprojectname -n ngspipe-rnaseq-basic
    ngspipedb runpipe myprojectname ngspipe-rnaseq-basic --report -db
    ngspipedb rundb server -url 0.0.0.0:8000
    
    A more detailed tutorial of how to use this toolkit can be found here: https://xuanblo.github.io/NGSPipeDb/
    """
    if args:
        click.echo(ctx.get_help())

def get_last_version_from_pypi(package, url_pattern):
    """Return version of package on pypi.python.org using json."""
    req = requests.get(url_pattern.format(package=package))
    version = parse('0')
    if req.status_code == requests.codes.ok:
        j = json.loads(req.text.encode(req.encoding))
        releases = j.get('releases', [])
        for release in releases:
            ver = parse(release)
            if not ver.is_prerelease:
                version = max(version, ver)
    if 'a' in str(version):
        dot_numbers_version = str(version).split('a')[0] # test alf version
    elif 'b' in str(version):
        dot_numbers_version = str(version).split('b')[0] # test beta version
    else:
        dot_numbers_version = str(version) # stable version 
    return dot_numbers_version

def new_version_check():
    URL_PATTERN = 'https://test.pypi.org/pypi/{package}/json'
    last_version = get_last_version_from_pypi('NGSPipeDb', URL_PATTERN)
    update_message = '''
WARNING: You are using NGSPipeDb version {current_version}; however, version {last_version} is available.
You should consider upgrading via the 'pip install -i https://test.pypi.org/simple/ ngspipedb={last_version}' command.
Please see Changelog for the latest changes: https://xuanblo.github.io/NGSPipeDb/changelog/
    '''.format(current_version=current_version, last_version=last_version)
    last_version_message = '''
You are using the lastest version {}
    '''.format(last_version)
    versions = Comparison(current_version, last_version)
    if versions.get_greater() == current_version or versions.get_greater() == None:
        print(last_version_message)
    else:
        print(update_message)

def ngspipedb_cli_main():
    cli.add_command(configure_env_group(args=None if sys.argv[3:] else ['-h']))
    cli.add_command(configure_download_group(args=None if sys.argv[2:] else ['-h']))
    cli.add_command(configure_ngspipe_group(args=None if sys.argv[3:] else ['-h']))
    cli.add_command(configure_ngsdb_group(args=None if sys.argv[3:] else ['-h']))
    cli.add_command(configure_startproject_group(args=None if sys.argv[3:] else ['-h']))
    #cli.add_command(configure_runserver_group(args=None if sys.argv[3:] else ['-h']))
    #cli.add_command(configure_info_group(args=None if sys.argv[3:] else ['-h']))
    new_version_check()
    cli(obj={}, args=None if sys.argv[1:] else ['-h'])

if __name__ == '__main__':
    ngspipedb_cli_main()