#!/usr/bin/env python3
# -*- coding:utf-8 -*-

from ruamel.yaml import YAML
import pathlib
import sys

def modify_config_yaml(args, ngspipedb_configfile, new_configfile):
    '''
    read yaml file, change a item, return str
    TODO: use abs path? or relative path
    '''
    modify_dict = {}
    if 'genomeFasta_path' in args.keys():
        modify_dict['genomeFasta_path'] = args['genomeFasta_path']
    if 'genomeAnno_path' in args.keys():
        modify_dict['genomeAnno_path'] = args['genomeAnno_path']
    if 'sample_path' in args.keys():
        modify_dict['sample_path'] = args['sample_path']
    if 'condition_path' in args.keys():
        modify_dict['condition_path'] = args['condition_path']
    if 'results_name' in args.keys():
        modify_dict['results_name'] = args['results_name']
    if 'rawreads_dir' in args.keys():
        modify_dict['rawreads_dir'] = args['rawreads_dir']
    if 'email_addr' in args.keys():
        modify_dict['email_addr'] = args['email_addr']
    if 'exp_path' in args.keys():
        modify_dict['exp_path'] = args['exp_path']
    if 'reads_prefix' in args.keys():
        modify_dict['read1Suffix'] = args['reads_prefix'].format('1')
        modify_dict['read2Suffix'] = args['reads_prefix'].format('2')
    #print(modify_dict)
    yaml = YAML()
    ngspipedb_configfile_path = pathlib.Path(ngspipedb_configfile)
    doc = yaml.load(ngspipedb_configfile_path)
    
    for k,v in modify_dict.items():
        doc[k] = v

    new_configfile_path = pathlib.Path(new_configfile)
    yaml.dump(doc, new_configfile_path)

    return 'done'

if __name__ == '__main__':
    args = {}
    ngspipedb_configfile = sys.argv[1]
    new_configfile = sys.argv[2]
    for i in sys.argv[3:]:
        k,v = i.split('=')
        args[k] = v
    #print(args)
    modify_config_yaml(args, ngspipedb_configfile, new_configfile) 
