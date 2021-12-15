import click
import time
import sys
import os

module_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(module_path)
from ngspipedbcli.download import download_main
from ngspipedbcli.runpipe import run_ngsdb, run_ngspipe, run_ngspipedb_one_step, run_server
from ngspipedbcli.project import start_project_main
from ngspipedbcli.env import env_create, env_list, env_pack, env_remove, env_unpack, env_update
from ngspipedbcli.common import *

@click.group()
@click.pass_context
def cli(ctx, args=None):
    """
    ngspipedb is a snakemake-based tool for reproducible next generation sequencing (NGS) data analysis and interactive web application auto-build.
    """
    pass

# /Users/zhangxuan/Work/Current_work2020-6-21/1.ngspipedb/packages/testdata_ngspipe-rnaseq-basic/

@cli.command()
@click.pass_context
@click.option('--projectname', default='xxx_analysis', help='[Required] (str) name your project', show_default=True, required=True)
@click.option('--pipename', help='[Required] pipeline', type=click.Choice(['ngspipe-rnaseq-basic']), required=True, show_default=True)
@click.option('--directory', help="[Required] (directory path) project directory, default output files under ~/ngspipedb_result/xxx_analysis", default='~/ngspipedb_result', show_default=True, type=click.STRING, required=True)
@click.option('--jobs', type=int, help="how many cpu to use (int)", default=1, show_default=True, required=True)
@click.option('--rawreadsdir', help="[Required] (directory path) raw reads directory", default='rawdata', show_default=True, type=click.STRING, required=True)
@click.option('--genomeFasta', help="[Required] (fasta file path) genome sequence file", default='genome/chr19.fa', show_default=True, type = click.STRING, required=True)
@click.option('--genomeAnno', help="[Required] (gff/gtf file path) genome annotation file", default='genome/GRCm38.83.chr19.gtf', show_default=True, type = click.STRING, required=True)
@click.option('--samplefile', help="[Required] (file path) samplefile", default='rawdata/sample.csv', show_default=True, type = click.STRING, required=True)
@click.option('--conditionfile', help="[Required] (file path) conditionfile", default='rawdata/condition.csv', show_default=True, type = click.STRING, required=True)
@click.option('--eggnogdir', help="[Required] (directory path) eggnog database directory", default='database/eggnog', show_default=True, type=click.STRING)
@click.option('--ontologyfile', help="[Required] (obo file path) (file path) gene ontology database file", default='database/gene_ontology/go-basic.obo', show_default=True, type=click.STRING)
@click.option('--emailaddr', help="[optional] example: xxx@qq.com", default='', type=click.STRING, show_default=True, required=False)
@click.option('--readsprefix', help="[Required] (str) reads prefix", default='_R{}.fq.gz', type=click.STRING, show_default=True, required=True)
@click.option('--resultdirname', help="[Required] (str) result directory name under proect name directory", default='results', type=click.STRING, show_default=True, required=True)
@click.option('--snaketype', help="`p`: print snakemake shell commands. `np`: Enable the dry run.", type=click.Choice(['np', 'p']), default='p')
@click.option('--report', help="[optional] generate html report", is_flag=True)
@click.option('--database', help="[optional] generate database", is_flag=True)
@click.option('--updateenv', help="[optional] install and update env auto", is_flag=True)
@click.option('--target', help="[optional] choose where to stop your pipeline", type=click.Choice(pipes_dict['ngspipe-rnaseq-basic']['steps'].keys()), default='all', show_default=True, required=True)
@click.option('--configfile', help="config file path", type=click.Path())
@click.option('--otherparams', help="[optional] other snakemake params Example: -f", default='')
@click.option('--printshell', is_flag=True, help="[optional] print shell commands, not real run")
def reference_based_rnaseq_analysis(ctx, **kargs):
    '''
    [very important]: users must make sure:
        1. path in Genomefasta/Genomeanno/Samplefile/Conditionfile inputbox exists
        2. [required] means inputbox must have values
        3. [optional] means inputbox can be empty
        4. `_R{}.fq.gz` in Reads prefix means ngspipedb will find *_R1.fq.gz and *_R2.fq.gz in Rawreadsdir
    '''
    click.echo(ctx.params)
    run_ngspipe(ctx.params)

@cli.command()
@click.pass_context
@click.option('--projectname', default='xxx_analysis', help='[Required] (str) name your project', show_default=True, required=True)
@click.option('--pipename', help='[Required] pipeline', type=click.Choice(['ngspipe-rnaseq-trinity']), required=True, show_default=True)
@click.option('--directory', help="[Required] (directory path) project directory, default output files under ~/ngspipedb_result/xxx_analysis", default='~/ngspipedb_result', show_default=True, type=click.STRING, required=True)
@click.option('--jobs', type=int, help="how many cpu to use (int)", default=1, show_default=True, required=True)
@click.option('--rawreadsdir', help="[Required] (directory path) raw reads directory", default='rawdata', show_default=True, type=click.STRING, required=True)
@click.option('--samplefile', help="[Required] (file path) samplefile", default='rawdata/sample.csv', show_default=True, type = click.STRING, required=True)
@click.option('--conditionfile', help="[Required] (file path) conditionfile", default='rawdata/condition.csv', show_default=True, type = click.STRING, required=True)
@click.option('--eggnogdir', help="[Required] (directory path) eggnog database directory", default='database/eggnog', show_default=True, type=click.STRING)
@click.option('--ontologyfile', help="[Required] (obo file path) (file path) gene ontology database file", default='database/gene_ontology/go-basic.obo', show_default=True, type=click.STRING)
@click.option('--emailaddr', help="[optional] example: xxx@qq.com", default='', type=click.STRING, show_default=True, required=False)
@click.option('--readsprefix', help="[Required] (str) reads prefix", default='_R{}.fq.gz', type=click.STRING, show_default=True, required=True)
@click.option('--resultdirname', help="[Required] (str) result directory name under proect name directory", default='results', type=click.STRING, show_default=True, required=True)
@click.option('--snaketype', help="`p`: print snakemake shell commands. `np`: Enable the dry run.", type=click.Choice(['np', 'p']), default='p')
@click.option('--maxmemeory', help="max memeory of trinity can use", type=click.STRING, default='2G', show_default=True, required=True)
@click.option('--report', help="[optional] generate html report", is_flag=True)
@click.option('--updateenv', help="[optional] install and update env auto", is_flag=True)
@click.option('--target', help="[optional] choose where to stop your pipeline", type=click.Choice(pipes_dict['ngspipe-rnaseq-trinity']['steps'].keys()), default='all', show_default=True, required=True)
@click.option('--configfile', help="config file path", type=click.Path())
@click.option('--otherparams', help="[optional] other snakemake params Example: -f", default='')
@click.option('--printshell', is_flag=True, help="[optional] print shell commands, not real run")
def denovo_rnaseq_analysis(ctx, **kargs):
    '''
    [very important]: users must make sure:
        1. path in Samplefile/Conditionfile inputbox exists
        2. [required] means inputbox must have values
        3. [optional] means inputbox can be empty
        4. `_R{}.fq.gz` in Reads prefix means ngspipedb will find *_R1.fq.gz and *_R2.fq.gz in Rawreadsdir
    '''
    click.echo(ctx.params)
    run_ngspipe(ctx.params)

@cli.command()
@click.pass_context
@click.option('--projectname', default='xxx_analysis', help='[Required] (str) name your project', show_default=True, required=True)
@click.option('--pipename', help='[Required] pipeline', type=click.Choice(['ngspipe-lncRNA']), required=True, show_default=True)
def lncRNA_analysis(ctx, **kargs):
    click.echo(ctx.params)
    run_ngspipe(ctx.params)

@cli.command()
@click.pass_context
@click.option('--projectname', default='xxx_analysis', help='[Required] (str) name your project', show_default=True, required=True)
@click.option('--pipename', help='[Required] pipeline', type=click.Choice(['ngspipe-chipseq']), required=True, show_default=True)
def chipseq_analysis(ctx, **kargs):
    click.echo(ctx.params)
    run_ngspipe(ctx.params)

@cli.command()
@click.pass_context
@click.option('--projectname', default='xxx_analysis', help='[Required] (str) name your project', show_default=True, required=True)
@click.option('--pipename', help='[Required] pipeline', type=click.Choice(['ngspipe-tnt']), required=True, show_default=True)
def tnt_insertion_analysis(ctx, **kargs):
    click.echo(ctx.params)
    run_ngspipe(ctx.params)

@cli.command()
@click.pass_context
@click.option('--projectname', default='xxx_analysis', help='[Required] (str) name your project', show_default=True, required=True)
@click.option('--pipename', help='[Required] pipeline', type=click.Choice(['ngspipe-resequencing']), required=True, show_default=True)
def resequcing_analysis(ctx, **kargs):
    click.echo(ctx.params)
    run_ngspipe(ctx.params)

if __name__ == '__main__':
    cli()