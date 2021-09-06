#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter, RawTextHelpFormatter
from textwrap import dedent
from rich import print
from rich.console import Console
from rich.table import Column, Table

from ngspipedbcli.download import download_main, print_all_avaiable_data
from ngspipedbcli.runpipe import run_ngsdb, run_ngspipe, run_test
from ngspipedbcli.project import startproject, copyproject
from ngspipedbcli.man import doc_download
from ngspipedbcli.server import run_server
from ngspipedbcli.env import env_main
from ngspipedbcli.common import __version__
from ngspipedbcli.common import *


COLORS = ['magenta', 'cyan']
console = Console()

class RichedArgParser(ArgumentParser):
    def _print_message(self, message, file=None):
        if message:
            if file is None:
                file = sys.stderr
            #file.write(message)
            from rich import print
            print(message)

def generate_parser():
    '''
    main ngspipedb_argparse
    Multi-level arg parse
    '''
    descr = dedent('''
        NGSPipeDb is an automated pipeline for parallel processing of huge next generation sequencing (NGS) data and database generation using Snakemake workflow which allows for ease of use, optimal speed, and a highly modular code that can be further added onto and customized by experienced users. It can be further divided into NGSPipe and NGSDb for individual usage.
    '''
    )
    additional_descr = dedent('''
        A more detailed tutorial of how to use this toolkit can be found here: https://xuanblo.github.io/NGSPipeDb/
    '''
    )
    p = RichedArgParser(
        prog='ngspipedb', 
        description= descr, 
        epilog=additional_descr, 
        formatter_class=RawDescriptionHelpFormatter,
        add_help=True,
    )
    p.add_argument("-V", "--version", action='version', version='ngspipedb {}'.format(__version__), help="Show the ngspipedb version number and exit.")
    subparsers = p.add_subparsers(
        title='[red]subcommands[/red]', 
        #metavar='command', 
        dest='cmd'
    )
    configure_parse_startproject(subparsers)
    configure_parse_copyproject(subparsers)
    configure_parse_downloaddb(subparsers)
    configure_parse_runngspipe(subparsers)
    configure_parse_runtest(subparsers)
    configure_parse_createcondaenv(subparsers)
    configure_parse_man(subparsers)
    configure_parse_check(subparsers)
    configure_parse_serve(subparsers)

    return p

# #############################################################################################
#
# sub-parsers
#
# #############################################################################################

def configure_parse_startproject(subparsers):
    '''
    startproject
    '''
    descr = dedent('''
        [green][u]startproject[/u] try to create a new directory for you, if the directory exists, it will not be created. 
        A typical NGS analysis project should include directories: Rawdata, Database, Script, Result.[/green]
    '''
    )
    additional_descr = dedent('''
        Example: 
        ngspipedb startproject -d run_snake_test -m rnaseq
    '''
    )
    p = subparsers.add_parser(
        'startproject', 
        help=descr,
        formatter_class=RawDescriptionHelpFormatter,
        epilog=additional_descr,
    )
    p.add_argument('-d', '--directory', help='workding directory', type=str, default='./', required=False)
    p.add_argument('-m', '--module', help='analysis module', type=str, default='rnaseq', required=True)
    p.add_argument('-n', '--dryrun', help='dry run', action='store_false')
    p.set_defaults(func=startproject)

def configure_parse_copyproject(subparsers):
    '''
    copyproject
    '''
    descr = dedent('''
        [blue]create a ngspipedb project auto[/blue]
    '''
    )
    additional_descr = dedent('''
        Example: 
        TODO
    '''
    )
    p = subparsers.add_parser(
        'copyproject', 
        help=descr,
        formatter_class=RawDescriptionHelpFormatter,
    )
    p.add_argument('--input', help='input project directory', type=str, required=True)
    p.add_argument('--output', help='output project directory', type=str, required=True)
    p.set_defaults(func=copyproject)

def configure_parse_downloaddb(subparsers):
    '''
    downloaddb
    '''
    descr = dedent('''
        [green]In order to annotate genes, including GO and KEGG, you need to download the corresponding database first.[/green]
    '''
    )
    additional_descr = dedent('''
        Example: 
        ngspipedb downloaddb -l http://www.liu-lab.com/ngspipedb/rnaseq_testdata.tar.gz -d run_snake_test
    '''
    )
    p = subparsers.add_parser(
        'downloaddb', 
        help= descr,
        formatter_class=RawDescriptionHelpFormatter,
        epilog=additional_descr,
    )
    p.add_argument('-l', '--url', help='url', required=True)
    p.add_argument('-d', '--directory', help='[green]download to directory[/green] (default: ./)', required=True, default='./')
    p.add_argument('-n', '--dryrun', help='dry run', action='store_false')
    p.set_defaults(func=download_main)

def configure_parse_runngspipe(subparsers):
    '''
    runtest
    '''
    descr = dedent('''
        [blue]pipeline can be rnaseq-basic, chipseq, resequencing[/blue]
    '''
    )
    additional_descr = dedent('''
        Example:
        TODO:
    '''
    )
    p = subparsers.add_parser(
        'runngspipe', 
        help=descr,
        formatter_class=RawDescriptionHelpFormatter,
        epilog=additional_descr, 
    )
    p.add_argument('-m', '--module', 
        help=dedent(
            '''
            [blue]pipeline can be rnaseq-basic, chipseq, resequencing[/blue]
            '''),
        choices=modules_dict.keys(),
        required=True,
        metavar='pipeline',
    )
    p.add_argument('-n', '--dryrun', help='dry run', action='store_false')
    p.set_defaults(func=run_ngspipe)

def configure_parse_runtest(subparsers):
    '''
    runtest
    '''
    descr = dedent('''
        [blue]test ngspipedb module[/blue]
    '''
    )
    additional_descr = dedent('''
        Example:
        ngspipedb runtest -m ngspipe-rnaseq-basic -d run_snake_test -t run_snake_test/testdata -j 1 -r
        ngspipedb runtest -m ngsdb -d run_snake_test -t run_snake_test/testdata -j 1
    '''
    )
    p = subparsers.add_parser(
        'runtest', 
        help=descr,
        formatter_class=RawDescriptionHelpFormatter,
        epilog=additional_descr, 
    )
    p.add_argument('-m', '--module', help='run pipeline with testdata', choices=['ngspipe-rnaseq-basic', 'ngsdb'])
    p.add_argument('-d', '--directory', help='working directory', required=True)
    p.add_argument('-t', '--testdata', help='testdata url', required=True)
    p.add_argument('-j', '--jobs', help='how many cpu to use', required=True, default=1, type=int)
    p.add_argument('-e', '--email', help='email address to recive', default='nobody', type=str)
    p.add_argument('-r', '--report', help='generate report', action='store_true')
    p.add_argument('-n', '--dryrun', help='dry run', action='store_false')
    p.set_defaults(func=run_test)

def configure_parse_createcondaenv(subparsers):
    '''
    createcondaenv
    '''
    descr = dedent('''
        [green]create ngspipe/ngsdb environment for ngs analysis & visualize[/green]
    '''
    )
    additional_descr = dedent('''
        Example:
        ngspipedb createcondaenv -m ngspipe-rnaseq-basic
        ngspipedb createcondaenv -m ngsdb
    '''
    )
    p = subparsers.add_parser(
        'createcondaenv', 
        help=descr,
        formatter_class=RawDescriptionHelpFormatter,
        epilog=additional_descr,
        )
    p.add_argument('-m', '--module', 
        help='create conda environment for module analysis',
        choices=modules_dict.keys(),
        required=True,
    )
    p.add_argument('-n', '--dryrun', help='dry run', action='store_false')
    p.set_defaults(func=env_main)

def configure_parse_man(subparsers):
    '''
    check
    '''
    descr = dedent('''
        check env, database, module
    '''
    )
    additional_descr = dedent('''
        [blue]read document[/blue]
    '''
    )
    p = subparsers.add_parser(
        'man', 
        help=additional_descr,
        formatter_class=RawDescriptionHelpFormatter,
    )
    p.add_argument('--ngspipe-rnaseq', help='[u cyan]read document for RNA-Seq analysis[/u cyan]', required=False)
    p.add_argument('--ngspipe-chipseq', help='read document  for Chip-Seq analysis', required=False)
    p.add_argument('--ngspipe-resequencing', help='read document  for Resequencing-Seq analysis', required=False)
    p.add_argument('--ngsdb', help='read document for database/website generate', required=False)
    p.set_defaults(func=doc_download)

def configure_parse_check(subparsers):
    '''
    check
    '''
    descr = dedent('''
        check env, database, module
    '''
    )
    additional_descr = dedent('''
        Example:
        ngspipedb check -h
    '''
    )
    p = subparsers.add_parser(
        'check', 
        help=descr, 
        epilog=additional_descr,
        formatter_class=RawDescriptionHelpFormatter,
    )
    p.add_argument(
        '--env',
    )

def configure_parse_serve(subparsers):
    '''
    serve
    '''
    descr = dedent('''
        run server
    '''
    )
    additional_descr = dedent('''
        Example:
        ngspipedb serve -d run_snake_test
    '''
    )
    p = subparsers.add_parser(
        'serve', 
        help=descr,
        formatter_class=RawDescriptionHelpFormatter,
        epilog=additional_descr,
    )
    p.add_argument('-d', '--directory', help='ngsdb result', required=True)
    p.add_argument('-n', '--dryrun', help='dry run', action='store_false')
    p.set_defaults(func=run_server)

def ngspipedb_cli_main():
    p = generate_parser()
    args = p.parse_args(None if sys.argv[1:] else ['-h'])
    args.func(args)
    #if args.verbosity:
    #    print('verbosity turned on')

if __name__ == '__main__':
    sys.exit(ngspipedb_cli_main())