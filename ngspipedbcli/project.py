#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from ngspipedbcli.common import *
from os.path import abspath, join

console = Console()

def start_ngspipe_project(args):
    '''
    TODO: add more pipeline
    TODO: add a flag file under .ngspipedb
    '''

    if args['directory']:
        working_dir = abspath(join(args['directory'], args['projectname']))
    else:
        working_dir = abspath(join('', args['projectname']))

    ngspipe_configfile = pipes_dict[args['pipeline']]['configfile']
    ngsdb_configfile = pipes_dict['ngsdb']['configfile']

    new_ngspipe_configfile = '{directory}/ngspipe_config.yaml'.format(directory = working_dir)
    new_ngsdb_configfile = '{directory}/ngsdb_config.yaml'.format(directory = working_dir)

    # directory structure
    create_directory_base_command = 'mkdir -p {directory}'.format(directory = working_dir)
    create_directory_rawdata_command = 'mkdir -p {directory}/rawdata'.format(directory = working_dir)
    touch_samplefile_command = 'touch {directory}/rawdata/sample.csv'.format(directory = working_dir)
    if args['pipeline'] == 'ngspipe-rnaseq-basic':
        touch_conditionfile_command = 'echo "sample_id,Sample,Tissue" >{directory}/rawdata/condition.csv'.format(directory = working_dir)
    create_directory_database_command = 'mkdir -p {directory}/database'.format(directory = working_dir)
    create_directory_genome_command = 'mkdir -p {directory}/genome'.format(directory = working_dir)
    create_directory_result_command = 'mkdir -p {directory}/results'.format(directory = working_dir)
    create_directory_flag_command = 'mkdir -p {directory}/.ngspipedb'.format(directory = working_dir)

    # configfile
    copy_ngspipe_configfile_command = 'cp {ngspipedb_configfile} {new_configfile}'.format(ngspipedb_configfile=ngspipe_configfile, new_configfile = new_ngspipe_configfile)
    copy_ngsdb_configfile_command = 'cp {ngspipedb_configfile} {new_configfile}'.format(ngspipedb_configfile=ngsdb_configfile, new_configfile = new_ngsdb_configfile)
    
    if args['printshell']:
        ngspipedb_print_command('create base directory', create_directory_base_command)
        ngspipedb_print_command('create rawdata directory under base', create_directory_rawdata_command)
        ngspipedb_print_command('touch samplefile under rawdata', touch_samplefile_command)
        if args['pipeline'] == 'ngspipe-rnaseq-basic':
            ngspipedb_print_command('touch conditionfile under rawdata', touch_conditionfile_command)
        ngspipedb_print_command('create database under base', create_directory_database_command)
        #ngspipedb_print_command('create result directory under base', create_directory_result_command)
        ngspipedb_print_command('copy ngspipe confilefile under base', copy_ngspipe_configfile_command)
        ngspipedb_print_command('copy ngsdb confilefile under base', copy_ngsdb_configfile_command)
        ngspipedb_print_command('create flag directory under base', create_directory_flag_command)
    else:
        #call_status_start_project_command = subprocess.call(start_project_command, shell=True, encoding='utf-8')
        call_status_create_directory_base_command = subprocess.call(create_directory_base_command, shell=True, encoding='utf-8')
        call_status_create_directory_rawdata_command = subprocess.call(create_directory_rawdata_command, shell=True, encoding='utf-8')
        call_status_touch_samplefile_command = subprocess.call(touch_samplefile_command, shell=True, encoding='utf-8')
        if args['pipeline'] == 'ngspipe-rnaseq-basic':
            call_status_touch_conditionfile_command = subprocess.call(touch_conditionfile_command, shell=True, encoding='utf-8')
        call_status_create_directory_database_command = subprocess.call(create_directory_database_command, shell=True, encoding='utf-8')
        call_status_create_directory_genome_command = subprocess.call(create_directory_genome_command, shell=True, encoding='utf-8')
        #call_status_create_directory_result_command = subprocess.call(create_directory_result_command, shell=True, encoding='utf-8')
        call_status_copy_configfile_command = subprocess.call(copy_ngspipe_configfile_command, shell=True, encoding='utf-8')
        call_status_copy_configfile_command = subprocess.call(copy_ngsdb_configfile_command, shell=True, encoding='utf-8')
        call_status_create_directory_flag_command = subprocess.call(create_directory_flag_command, shell=True, encoding='utf-8')

def start_ngsdb_project(args):
    '''
    TODO: add more pipeline
    TODO: add a flag file under .ngspipedb
    '''

    if args['directory']:
        working_dir = abspath(join(args['directory'], args['projectname']))
    else:
        working_dir = abspath(join('', args['projectname']))

    ngspipedb_configfile = pipes_dict[args['pipeline']]['configfile']
    new_configfile = '{directory}/ngsdb_config.yaml'.format(directory = working_dir)
    # directory structure
    create_directory_base_command = 'mkdir -p {directory}'.format(directory = working_dir)
    create_directory_rawdata_command = 'mkdir -p {directory}/rawdata'.format(directory = working_dir)
    touch_samplefile_command = 'touch {directory}/rawdata/sample.csv'.format(directory = working_dir)
    if args['pipeline'] == 'ngspipe-rnaseq-basic':
        touch_conditionfile_command = 'echo "sample_id,Sample,Tissue" >{directory}/rawdata/condition.csv'.format(directory = working_dir)
    create_directory_database_command = 'mkdir -p {directory}/database'.format(directory = working_dir)
    create_directory_genome_command = 'mkdir -p {directory}/genome'.format(directory = working_dir)
    create_directory_result_command = 'mkdir -p {directory}/results'.format(directory = working_dir)
    create_directory_flag_command = 'mkdir -p {directory}/.ngspipedb'.format(directory = working_dir)
    # configfile
    copy_configfile_command = 'cp {ngspipedb_configfile} {new_configfile}'.format(ngspipedb_configfile=ngspipedb_configfile, new_configfile = new_configfile)
    
    if args['printshell']:
        ngspipedb_print_command('create base directory', create_directory_base_command)
        ngspipedb_print_command('create rawdata directory under base', create_directory_rawdata_command)
        ngspipedb_print_command('touch samplefile under rawdata', touch_samplefile_command)
        if args['pipeline'] == 'ngspipe-rnaseq-basic':
            ngspipedb_print_command('touch conditionfile under rawdata', touch_conditionfile_command)
        ngspipedb_print_command('create database under base', create_directory_database_command)
        #ngspipedb_print_command('create result directory under base', create_directory_result_command)
        ngspipedb_print_command('copy confilefile under base', copy_configfile_command)
        ngspipedb_print_command('create flag directory under base', create_directory_flag_command)
    else:
        #call_status_start_project_command = subprocess.call(start_project_command, shell=True, encoding='utf-8')
        call_status_create_directory_base_command = subprocess.call(create_directory_base_command, shell=True, encoding='utf-8')
        call_status_create_directory_rawdata_command = subprocess.call(create_directory_rawdata_command, shell=True, encoding='utf-8')
        call_status_touch_samplefile_command = subprocess.call(touch_samplefile_command, shell=True, encoding='utf-8')
        if args['pipeline'] == 'ngspipe-rnaseq-basic':
            call_status_touch_conditionfile_command = subprocess.call(touch_conditionfile_command, shell=True, encoding='utf-8')
        call_status_create_directory_database_command = subprocess.call(create_directory_database_command, shell=True, encoding='utf-8')
        call_status_create_directory_genome_command = subprocess.call(create_directory_genome_command, shell=True, encoding='utf-8')
        #call_status_create_directory_result_command = subprocess.call(create_directory_result_command, shell=True, encoding='utf-8')
        call_status_copy_configfile_command = subprocess.call(copy_configfile_command, shell=True, encoding='utf-8')
        call_status_create_directory_flag_command = subprocess.call(create_directory_flag_command, shell=True, encoding='utf-8')

def start_project_main(args):
    if args['pipeline'] == 'ngsdb':
        start_ngsdb_project(args)
    else:
        start_ngspipe_project(args)

if __name__ == '__main__':
    print(__name__)