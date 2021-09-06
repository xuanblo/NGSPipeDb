#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import time
from rich import print
from rich.console import Console
from rich.progress import track
from os.path import join

from ngspipedbcli.common import *

console = Console()

def get_remote_file_size(url):
    import urllib
    import requests
    download = urllib.request.urlopen(url)
    return download.headers['content-length']

def download_with_progress_bar(url, target_file_name):
    from urllib.request import urlopen
    response = urlopen(url)
    chunk_file_size = 16 * 1024
    remote_file_size = int(get_remote_file_size(url))
    fetch_file_size = 0
    with open(target_file_name, 'wb') as f:
        #while fetch_file_size < int(remote_file_size):
        for step in track(range(round(remote_file_size/chunk_file_size)+1), description="Downloading {}".format(target_file_name)):
            chunk = response.read(chunk_file_size)
            if not chunk:
                break
            fetch_file_size += chunk_file_size
            f.write(chunk)

def test_paramiko():
    pass

def cloud_storage():
    pass

def _wget():
    'wget -O {output.ref_gff} https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore\&report=gff3\&id=$accession --tries 10 | tee {log} 2> {log}'

def download_main(args):
    '''
    download
    '''
    if args['list']:
        # print a table with rich
        print_rich_markdown_data_table()
        return
    pipeline_dict = call_snakemake_module(args['pipeline'])

    # outputfile
    datatype = args['datatype']
    pipeline = args['pipeline']
    platform = args['platform']
    target_file_name = []
    if datatype == 'env' and platform=='linux':
        target_file_name.append(pipeline_dict['packenv'][0])
    elif datatype == 'env' and platform=='osx':
        target_file_name.append(pipeline_dict['packenv'][1])
    elif datatype == 'database':
        target_file_name = pipeline_dict['database']
    elif datatype == 'testdata':
        target_file_name.append(pipeline_dict['testdata'])
    comands = []
    for i in target_file_name:
        command_download_data = "wget -P {directory} {url}".format(directory=args['directory'], url=join(dataurl, i))
        comands.append(command_download_data)
    
    if args['printshell']:
        for i in comands:
            ngspipedb_print_command('downloading...', i)
    else:
        for i in comands:
            call_status_command_download_testdata = subprocess.call(i, shell=True, encoding='utf-8')

if __name__ == '__main__':
    args = {}
    download_main(args)