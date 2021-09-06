#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from ngspipedbcli.common import *

def run_server(args):
    '''
    run ngsdb result
    '''

    # step1, create a conda env by gived name
    # print
    run_ngsdb_django_database_command = 'source activate ngsdb && cd {directory} && python results/ngsdb/manage.py runserver && conda deactivate'.format(directory=args.directory)
    ngspipedb_print_command(run_ngsdb_django_database_command)

    # step2, run
    if args.dryrun:
        p = subprocess.Popen(run_ngsdb_django_database_command, shell=True, encoding='utf-8')
        try:
            p.wait()
        except KeyboardInterrupt:
            try:
                p.terminate()
            except OSError:
                pass
            p.wait()