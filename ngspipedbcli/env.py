#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from ngspipedbcli.common import *
from os.path import join

"""
def env_ngspipe_rnaseq_basic(args):
    '''
    this is a two step create conda env, replaced by one step
    '''

    # step1, create a conda env by gived name
    # print
    create_conda_env_command = 'mamba create -c conda-forge -c bioconda -y --name {name} snakemake=5.30.2 python=3.8 seqkit=0.14.0'.format(name = args['name'])
    
    # run
    if args['printshell']:
        ngspipedb_print_command('create environment', create_conda_env_command)
    else:
        call_status_create_conda_env_command = subprocess.call(create_conda_env_command, shell=True, encoding='utf-8')

    # step2, update bioinformatics tools
    # print
    modules_dict = call_snakemake_module(args['name'])
    update_conda_env_command = 'mamba env update --name {name} --file {env} --prune'.format(name = args['name'], env=modules_dict['env_path'])
    
    # run
    if args['printshell']:
        ngspipedb_print_command('update environment', update_conda_env_command)
    else:
        call_status_update_conda_env_command = subprocess.call(update_conda_env_command, shell=True, encoding='utf-8')
"""

def env_create(args):
    '''
    create conda environment for pipeline
    '''
    pipe_dict = call_snakemake_module(args['name'])
    create_conda_env_command = 'mamba env create --name {name} --file {env_path}'.format(name = args['name'], env_path=pipe_dict['env_path'])

    if args['printshell']:
        ngspipedb_print_command('create environment', create_conda_env_command)
    else:
        call_status_update_conda_env_command = subprocess.call(create_conda_env_command, shell=True, encoding='utf-8')

def env_list(args):
    '''
    conda env list
    '''
    list_conda_env_command = 'conda env list'
    call_status_list_conda_env_command = subprocess.run(list_conda_env_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, encoding='utf-8', universal_newlines=True)
    env_dict = dict()
    for col, line in enumerate(call_status_list_conda_env_command.stdout.split('\n')):
        import re
        matObj = re.match('(\S+)[^/]+(/.*)', line)
        #print(col, line)
        #if matObj:
        #    print(matObj.groups())
        #else:
        #    print('------')
        if matObj:
            env_name = matObj.group(1)
            if env_name in ngspipedb_envs_set:
                line = line.replace(env_name, '[red]'+env_name+'[red]', 1)
        #ngspipedb_print_default_stdout(line)
        print(line)

def env_pack(args):
    if args['directory']:
        target_name = join(args['directory'], '{name}_{platform}.tar.gz'.format(name=args['name'], platform=args['platform']))
    else:
        target_name = '{name}_{platform}.tar.gz'.format(name=args['name'], platform=args['platform'])
    pack_conda_env_command = 'mamba pack -n {name} -o {target_name}.tar.gz && chmod +r {target_name}.tar.gz'.format(name=args['name'], target_name=target_name)
    '''TODO: chmod +r xxx.tar.gz'''
    if args['printshell']:
        ngspipedb_print_command('pack environment', pack_conda_env_command)
    else:
        call_status_pack_conda_env_command = subprocess.call(pack_conda_env_command, shell=True, encoding='utf-8')

def env_unpack(args):
    name = os.path.basename(args['file']).split('_')[0]
    unpack_conda_env_command = 'mkdir -p {condaenvdir}/{name} && tar -xzf {filename} -C {condaenvdir}/{name} && source {condaenvdir}/{name}/bin/activate && conda-unpack && source {condaenvdir}/{name}/bin/deactivate'.format(condaenvdir=args['condaenvdir'], name=name, filename=args['file'])
    if args['printshell']:
        ngspipedb_print_command('unpack environment', unpack_conda_env_command)
    else:
        call_status_unpack_conda_env_command = subprocess.call(unpack_conda_env_command, shell=True, encoding='utf-8')

def env_remove(args):
    if args['name']:
        remove_conda_env_command = 'conda env remove -n {name}'.format(name=args['name'])
        if args['printshell']:
            ngspipedb_print_command('remove conda name...', remove_conda_env_command)
        else:
            call_status_remove_conda_env_command = subprocess.call(remove_conda_env_command, shell=True, encoding='utf-8')
    if args['path']:
        remove_conda_env_command = 'conda env remove -p {path}'.format(name=args['path'])
        if args['printshell']:
            ngspipedb_print_command('remove conda path...', remove_conda_env_command)
        else:
            call_status_remove_conda_env_command = subprocess.call(remove_conda_env_command, shell=True, encoding='utf-8')

def env_update(args):
    if args['name']:
        update_conda_env_command = 'mamba env update -n {name} -f {yaml}'.format(name=args['name'], yaml=pipes_dict[args['pipename']]['env_path'])
        if args['printshell']:
            ngspipedb_print_command('update conda name...', update_conda_env_command)
        else:
            call_status_remove_conda_env_command = subprocess.call(update_conda_env_command, shell=True, encoding='utf-8')
    if args['path']:
        update_conda_env_command = 'mamba env update -p {path} -f {yaml}'.format(name=args['path'], yaml=pipes_dict[args['pipename']]['env_path'])
        if args['printshell']:
            ngspipedb_print_command('remove conda path...', update_conda_env_command)
        else:
            call_status_remove_conda_env_command = subprocess.call(update_conda_env_command, shell=True, encoding='utf-8')

if __name__ == '__main__':
    pass