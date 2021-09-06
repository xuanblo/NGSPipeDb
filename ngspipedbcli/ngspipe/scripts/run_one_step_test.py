import os
import logging

from optparse import OptionParser

usage = '''
    One time installation of components necessary for an individual use.

    The default result was stored  under current directory.

    if -e not 指定, env will not be installed. conda path or xxx path.

        The gene_matrix_file.

    Alternatively, you can specify the filenames directly with
    :option:`-p`/:option:`--condapath` and
    :option:`-n`/:option:`--condaname`
    :option:`-d`/:option:`--testdatapath`
    :option:`-t`/:option:`--pathform` option.

    Example::

        python run_one_step_test.py -p linux -d testdata -e condapath -g gene_matrix_file
'''

parser = OptionParser(usage)
parser.add_option('-t', '--platform', dest='type', help='mac or linux', metavar="STR", action = "store", type="string")
parser.add_option("-d","--testdatadir", dest="stringtiedir", help="download test file", metavar="DIR", action = "store", type="string")
parser.add_option("-n","--name", dest="condaname", help="download ngspipedb env", metavar="FILE", action = "store", type="string")
parser.add_option("-p","--path", dest="condapath", help="download ngspipedb env", metavar="FILE", action = "store", type="string")

(options, args)=parser.parse_args()


import subprocess


def install_sh():
    try:
        retcode = subprocess.call("pip install sh", shell=True)
        return retcode
    except OSError as e:
        return "Execution failed:", e


try:
    import sh
except ImportError:
    install_sh()
    import sh

# ps -auxc | grep nginx
def is_nginx_running():
    r = sh.grep(sh.ps("-auxc"), "nginx", _ok_code=[1, 2, 3])
    return r.exit_code == 0


def install_nginx():
    if not sh.which("nginx"):
        print "nginx not exist, will install"
        sh.apt_get("install", "nginx", "-y")
    else:
        print "nginx has installed"


def start_nginx():
    r = sh.service("nginx", "start", _ok_code=[1, 2, 3])
    if r.exit_code == 0:
        print "start success"
    else:
        print "start failed"


if __name__ == "__main__":
    if not is_nginx_running():
        install_nginx()
        start_nginx()
    else:
        print "nginx is running"
        
import coloredlogs, logging

# Create a logger object.
logger = logging.getLogger(__name__)

# By default the install() function installs a handler on the root logger,
# this means that log messages from your code and log messages from the
# libraries that you use will all show up on the terminal.
coloredlogs.install(level='DEBUG')

# If you don't want to see log messages from libraries, you can pass a
# specific logger object to the install() function. In this case only log
# messages originating from that logger will show up on the terminal.
coloredlogs.install(level='DEBUG', logger=logger)

logger.debug("this is a debugging message")
logger.info("this is an informational message")
logger.warning("this is a warning message")
logger.error("this is an error message")
logger.critical("this is a critical message")

import sh
aa=sh.ls()
#输出：(0, 'total 0\n-rw-rw-r-- 1 roaddb roaddb 0 Dec 11 10:09 a.txt\n-rw-rw-r-- 1 roaddb roaddb 0 Dec 11 10:09 b.txt') 
logger.info(aa)

#os.system("Miniconda3-latest-Linux-x86_64.sh && bash /tmp Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3")
#wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/Miniconda3-latest-Linux-x86_64.sh && bash /tmp/Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3