
import sys
import os
from os.path import join
from ngspipedbcli.common import *

def gui_main(args):
    #step1 = 'export FLASK_ENV=development'
    #step2 = 'export FLASK_APP={}'.format(join(guiwebapp_dir, 'runpipeapp.py'))
    #step2 = 'export FLASK_APP=ngspipedb:guiwebapp:runpipeapp'
    env_dict = {
        'FLASK_ENV': 'development',
        'FLASK_APP': join(guiwebapp_dir, 'runpipeapp.py')
    }
    step3 = 'flask run -h {} -p {}'.format(args['ip'], args['port'])
    if args['printshell']:
        #ngspipedb_print_command('flask env', step1)
        #ngspipedb_print_command('flask env', step2)
        ngspipedb_print_command('flask run', step3)
    else:
        print('start flask server')
        #call_step1 = subprocess.call(step1, shell=True, encoding='utf-8')
        #call_step2 = subprocess.call(step2, shell=True, encoding='utf-8')
        os.putenv('FLASK_ENV', 'development')
        os.putenv('FLASK_APP', join(guiwebapp_dir, 'runpipeapp.py'))
        call_step3 = subprocess.call(step3, shell=True, encoding='utf-8')