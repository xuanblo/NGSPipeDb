#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from os.path import abspath, join, exists
import sys
from ngspipedbcli.common import *

def make_parmas_flat(args):
    flat_params_list = []
    if 'genomefasta' in args.keys() and args['genomefasta']:
        flat_params_list.append('genomeFasta_path={}'.format(abspath(args['genomefasta'])))
    if 'genomeanno' in args.keys() and args['genomeanno']:
        flat_params_list.append('genomeAnno_path={}'.format(abspath(args['genomeanno'])))
    if 'samplefile' in args.keys() and args['samplefile']:
        flat_params_list.append('sample_path={}'.format(abspath(args['samplefile'])))
    if 'rawreadsdir' in args.keys() and args['rawreadsdir']:
        flat_params_list.append('rawreads_dir={}'.format(abspath(args['rawreadsdir'])))
    if 'conditionfile' in args.keys() and args['conditionfile']:
        flat_params_list.append('condition_path={}'.format(abspath(args['conditionfile'])))
    if 'email_addr' in args.keys() and args['email_addr']:
        flat_params_list.append('email_addr={}'.format(args['email_addr']))
    if 'resultdirname' in args.keys() and args['resultdirname']:
        flat_params_list.append('results_name={}'.format(args['resultdirname']))
    if 'exp_profile' in args.keys() and args['exp_profile']:
        flat_params_list.append('exp_path={}'.format(abspath(args['exp_profile'])))
    if 'reads_prefix' in args.keys() and args['reads_prefix']:
        flat_params_list.append('reads_prefix={}'.format(args['reads_prefix']))
    return ' '.join(flat_params_list)

def check_config_paths(workding_directory, configfile):
    '''
    all path in configfile exists and file_size>0
    '''
    from ruamel.yaml import YAML
    import pathlib
    yaml = YAML()
    ngspipedb_configfile_path = pathlib.Path(configfile)
    doc = yaml.load(ngspipedb_configfile_path)

    for k,v in doc.items():
        # check path exists and empty
        v_abs = join(workding_directory, v)
        if k.endswith('_path'):
            #print(k, v)
            # chech if v is a abs path
            if not v.startswith('/'):
                v = v_abs
                #print(v)
            if os.path.isfile(v) and os.path.getsize(v) > 0:
                pass
            else:
                sys.stderr.write('param: {k} is required, however its path: {v} is not exists or empty!\n'.format(k=k, v=v))
                sys.exit(-1)
        # check dir exists and empty
        if k.endswith('_dir'):
            #print(k, v)
            # chech if v is a abs path
            if not v.startswith('/'):
                v = v_abs
            if os.path.isdir(v):
                # if fastq file in this directory then create links
                if os.listdir(v):
                    pass
                else:
                    sys.stderr.write('param: {k} is required, however no fastq file in dir: {v} is not exists\n'.format(k=k, v=v))
            else:
                sys.stderr.write('param: {k} is required, however its dir: {v} is not exists\n'.format(k=k, v=v))
                sys.exit(-1)

def save_command(workshpath):
    '''
    TODO: if pipename is not given, we will guess by .ngspipedb
    '''
    pass

def check_ngspipedb_conda_env(conda_env):
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
            env_dict[matObj.group(1)] = matObj.group(2)
    if conda_env in env_dict.keys():
        print('environment exists.')
        return True
    else:
        print('environment not exists. trying to install automaticly.')
        return False

def run_server(args):
    '''
    run ngsdb result
    '''

    # step1, create a conda env by gived name
    # print
    run_ngsdb_django_database_command = 'source activate ngsdb && python {managepypath} runserver {urlport} && conda deactivate'.format(managepypath=args['managepy'], urlport=args['urlport'])
    if args['printshell']:
        ngspipedb_print_command('run server', run_ngsdb_django_database_command)
    else:
        # step2, run
        p = subprocess.Popen(run_ngsdb_django_database_command, shell=True, encoding='utf-8')
        try:
            p.wait()
        except KeyboardInterrupt:
            try:
                p.terminate()
            except OSError:
                pass
            p.wait()

def run_ngspipe(args):
    '''
    run test data 
    '''

    pipe_dict = call_snakemake_module(args['pipename'])
    conda_env = pipe_dict['env_name']

    project_name = args['projectname']
    workding_directory = abspath(join(args['directory'], project_name))

    run_info = '''
#-------------------------------------------#
pipe_name: {}
project_name: {}
workding_directory: {}
conda_env_name: {}
#-------------------------------------------#
    '''.format(args['pipename'], project_name, workding_directory, conda_env)
    ngspipedb_print_rich_stdout(run_info)

    # check environment
    if not check_ngspipedb_conda_env(conda_env):
        install_env_command = 'python -m ngspipedbcli env create -n {}'.format(conda_env)
        if args['printshell']:
            ngspipedb_print_command('install environment', install_env_command)
        else:
            run_status_install_env_command = subprocess.run(install_env_command, shell=True, encoding='utf-8')

    # create directory
    if not args['directory']:
        print('using current directory')
        create_project_command = 'python -m ngspipedbcli startproject {projectname} -n {pipe_name}'.format(projectname=args['projectname'], pipe_name=args['pipename'])
        if args['printshell']:
            ngspipedb_print_command('create a rnaseq-basic project directory structure', create_project_command)
        else:
            run_status_create_project_command = subprocess.run(create_project_command, shell=True, encoding='utf-8')
    elif exists(join(workding_directory, '.ngspipedb')):
        print('use a project directory pre-created')
    else:
        print('use a giving directory')
        create_project_command = 'python -m ngspipedbcli startproject {projectname} -n {pipe_name} -d {project_dir}'.format(projectname=args['projectname'], pipe_name=args['pipename'], project_dir=args['directory'])
        if args['printshell']:
            ngspipedb_print_command('create a rnaseq-basic project directory structure', create_project_command)
        else:
            run_status_create_project_command = subprocess.run(create_project_command, shell=True, encoding='utf-8')

    # modify configfile
    ngspipedb_configfile = '{working_dir}/ngspipe_config.yaml'.format(working_dir=workding_directory)
    if not args['configfile']:
        print('using basic configfile')
        configfile = ngspipedb_configfile
    else:
        print('using custom configfile')
        configfile = abspath(args['configfile'])

    if args['genomefasta'] or args['genomeanno'] or args['email_addr'] or args['reads_prefix'] or args['resultdirname']:
        from tempfile import NamedTemporaryFile
        tmp_configfile = NamedTemporaryFile(mode="w+", delete=True)
        tmp_configfile.write('for temp configfile')
        tmp_configfile.close()
        #print(tmp_configfile.name)
        script_path = '{}/scripts/modify_config_yaml.py'.format(ngspipe_dir)

        modify_configfile_command = 'python {script} {configfile} {tmp_configfile} {params_flat}'.format(script=script_path, configfile=configfile, tmp_configfile = tmp_configfile.name, params_flat=make_parmas_flat(args))

        move_comfigfile_command = 'mv {tmp_configfile} {ngspipedb_configfile}'.format(tmp_configfile=tmp_configfile.name, ngspipedb_configfile=ngspipedb_configfile)

        run_ngspipe_rnaseq_basic_command = 'source activate {env_name} && snakemake -s {snakefile} --configfile {config} --directory {working_dir} --rerun-incomplete --scheduler greedy --nolock --jobs {cores} -{snaketype} {otherparams} && conda deactivate'.format(env_name=args['pipename'], snakefile = pipe_dict['snakefile'], config=ngspipedb_configfile, working_dir=workding_directory, cores=args['jobs'], snaketype=args['snaketype'], otherparams=args['otherparams'])

        # step1 runpipe
        if args['printshell']:
            ngspipedb_print_command('modify configfile', modify_configfile_command)
            ngspipedb_print_command('update configfile', move_comfigfile_command)
            ngspipedb_print_command('run pipeline', run_ngspipe_rnaseq_basic_command)
        else:
            call_status_modify_configfile_command = subprocess.call(modify_configfile_command, shell=True, encoding='utf-8')
            call_status_move_comfigfile_command = subprocess.call(move_comfigfile_command, shell=True, encoding='utf-8')
            # check configfile
            check_config_paths(workding_directory, ngspipedb_configfile)
            print('run {pipeline} analysis'.format(pipeline=args['pipename']))
            call_status_activate_conda_env_command = subprocess.call(run_ngspipe_rnaseq_basic_command, shell=True, encoding='utf-8')
            

    # step2 ngspipe rnaseq basic report
    if args['report']:
        report_ngspipe_rnaseq_basic_command = 'source activate {env_name} && snakemake -s {snakefile} --configfile {configfile} --directory {working_dir} --report {report} -{snaketype} && conda deactivate'.format(env_name=args['pipename'], snakefile = pipe_dict['snakefile'], configfile=configfile, working_dir=workding_directory, report=join(workding_directory, args['resultdirname'], 'report', 'index.html'), snaketype=args['snaketype'])
        
        if args['printshell']:
            ngspipedb_print_command('generate report', report_ngspipe_rnaseq_basic_command)
        else:
            if exists(join(workding_directory, '.ngspipedb', 'run_all.runned')):
                print('run {pipeline} report'.format(pipeline=args['pipename']))
                call_status_report_ngspipe_rnaseq_basic_command = subprocess.call(report_ngspipe_rnaseq_basic_command, shell=True, encoding='utf-8')
            else:
                print('{} is not exists, may be pipeline haven\'t been run.'.format(join(workding_directory, '.ngspipedb', 'run_all.runned')))
    
    # step3 generate database
    if args['database']:
        if args['directory']:
            run_ngsdb_rnaseq_basic_command = 'python -m ngspipedbcli rundb build {projectname} -d {directory} --resultdirname {resultname} -exp {exp_path} --genomeFasta {genomeFasta} --genomeAnno {genomeAnno} --snaketype {snaketype}'.format(projectname=args['projectname'], directory=args['directory'], resultname=args['resultdirname'], exp_path=join(workding_directory, args['resultdirname'], 'ngspipe_result', 'quantify/quantify_by_stringtie/gene_fpkm_all_samples.tsv'), snaketype=args['snaketype'], genomeFasta=args['genomefasta'], genomeAnno=args['genomeanno'])
        else:
            run_ngsdb_rnaseq_basic_command = 'python -m ngspipedbcli rundb build {projectname} --resultdirname {resultname} -exp {exp_path} --genomeFasta {genomeFasta} --genomeAnno {genomeAnno} --snaketype {snaketype}'.format(projectname=args['projectname'], directory=args['directory'], resultname=args['resultdirname'], exp_path=join(workding_directory, args['resultdirname'], 'ngspipe_result', 'quantify/quantify_by_stringtie/gene_fpkm_all_samples.tsv'), snaketype=args['snaketype'], genomeFasta=args['genomefasta'], genomeAnno=args['genomeanno'])
        if args['printshell']:
            ngspipedb_print_command('build database', run_ngsdb_rnaseq_basic_command)
        else:
            call_status_run_ngsdb_rnaseq_basic_command = subprocess.call(run_ngsdb_rnaseq_basic_command, shell=True, encoding='utf-8')
            print('Vist database by command: ngspipedb rundb serve -m {managepypath} -up {urlport}'.format(managepypath=join(workding_directory, args['resultdirname'], 'ngsdb_code', 'manage.py'), urlport='0.0.0.0:8909'))

def run_ngsdb(args):
    '''
    generate database
    '''
    pipename = 'ngsdb'
    pipe_dict = call_snakemake_module(pipename)
    conda_env = pipe_dict['env_name']

    project_name = args['projectname']
    workding_directory = abspath(join(args['directory'], project_name))

    run_info = '''
#-------------------------------------------#
pipe_name: {}
project_name: {}
workding_directory: {}
conda_env_name: {}
demo_configfile: {}
#-------------------------------------------#
    '''.format(pipename, project_name, workding_directory, conda_env, abspath(pipe_dict['configfile']))
    ngspipedb_print_rich_stdout(run_info)

    # check environment
    if not check_ngspipedb_conda_env(conda_env):
        install_env_command = 'python -m ngspipedbcli env create -n {}'.format(conda_env)
        if args['printshell']:
            ngspipedb_print_command('install environment', install_env_command)
        else:
            run_status_install_env_command = subprocess.run(install_env_command, shell=True, encoding='utf-8')

    # create directory
    if not args['directory']:
        print('using current directory')
        create_project_command = 'python -m ngspipedbcli startproject {projectname} -n {pipe_name}'.format(projectname=args['projectname'], pipe_name=pipename)
        if args['printshell']:
            ngspipedb_print_command('create a rnaseq-basic project directory structure', create_project_command)
        else:
            run_status_create_project_command = subprocess.run(create_project_command, shell=True, encoding='utf-8')
    elif exists(join(workding_directory, '.ngspipedb')):
        print('use a project directory pre-created')
    else:
        print('use a giving directory')
        create_project_command = 'python -m ngspipedbcli startproject {projectname} -n {pipe_name} -d {project_dir}'.format(projectname=args['projectname'], pipe_name=pipename, project_dir=args['directory'])
        if args['printshell']:
            ngspipedb_print_command('create a rnaseq-basic project directory structure', create_project_command)
        else:
            run_status_create_project_command = subprocess.run(create_project_command, shell=True, encoding='utf-8')

    # modify configfile
    ngspipedb_configfile = '{working_dir}/ngsdb_config.yaml'.format(working_dir=workding_directory)
    
    if not args['configfile']:
        print('using basic configfile')
        if exists(ngspipedb_configfile):
            configfile = ngspipedb_configfile
        else:
            configfile = abspath(pipes_dict['ngsdb']['configfile'])
    else:
        print('using custom configfile')
        configfile = abspath(args['configfile'])

    if args['genomefasta'] or args['genomeanno'] or args['email_addr'] or args['exp_profile'] or args['resultdirname']:
        from tempfile import NamedTemporaryFile
        tmp_configfile = NamedTemporaryFile(mode="w+", delete=True)
        tmp_configfile.write('for temp configfile')
        tmp_configfile.close()
        #print(tmp_configfile.name)
        script_path = '{}/scripts/modify_config_yaml.py'.format(ngspipe_dir)

        modify_configfile_command = 'python {script} {configfile} {tmp_configfile} {params_flat}'.format(script=script_path, configfile=configfile, tmp_configfile = tmp_configfile.name, params_flat=make_parmas_flat(args))

        move_comfigfile_command = 'mv {tmp_configfile} {ngspipedb_configfile}'.format(tmp_configfile=tmp_configfile.name, ngspipedb_configfile=ngspipedb_configfile)

        run_ngsdb_command = 'source activate {env_name} && snakemake -s {snakefile} --configfile {config} --directory {working_dir} --rerun-incomplete --scheduler greedy --nolock --jobs {cores} -{snaketype} {otherparams} && conda deactivate'.format(env_name=pipename, snakefile = pipe_dict['snakefile'], config=ngspipedb_configfile, working_dir=workding_directory, cores=args['jobs'], snaketype=args['snaketype'], otherparams=args['otherparams'])

        # step1 runpipe
        if args['printshell']:
            ngspipedb_print_command('modify configfile', modify_configfile_command)
            ngspipedb_print_command('update configfile', move_comfigfile_command)
            ngspipedb_print_command('run pipeline', run_ngsdb_command)
        else:
            call_status_modify_configfile_command = subprocess.call(modify_configfile_command, shell=True, encoding='utf-8')
            call_status_move_comfigfile_command = subprocess.call(move_comfigfile_command, shell=True, encoding='utf-8')
            # check configfile
            check_config_paths(workding_directory, ngspipedb_configfile)
            print('run {pipeline} build'.format(pipeline=pipename))
            call_status_activate_conda_env_command = subprocess.call(run_ngsdb_command, shell=True, encoding='utf-8')
            

def run_ngspipedb_one_step(args):
    '''
    rnaseq, chipseq
    '''
    module_dict = call_snakemake_module(args.module)
    # step1 activate environment
    activate_conda_env_command = "conda activate {}".format(args.module)
    ngspipedb_print_command(activate_conda_env_command)
    #call_status_activate_conda_env_command = subprocess.call(activate_conda_env_command, shell=True, encoding='utf-8')
    # step2 run snakemake
    if args.dryrun:
        run_type = '-np'
    else:
        run_type = '-p'
    dryrun_ngspipe_command = 'snakemake -s {snakefile} --configfile {config} {run_type}'.format(snakefile = module_dict['snakefile'], config=module_dict['config'], run_type=run_type)
    ngspipedb_print_command(dryrun_ngspipe_command)
    call_status_activate_conda_env_command = subprocess.call(activate_conda_env_command, shell=True, encoding='utf-8')

    ngspipedb_configfile = pipes_dict['ngspipe-rnaseq-basic']['configfile']
    new_configfile = '{directory}/config.yaml'.format(directory = args['directory'])