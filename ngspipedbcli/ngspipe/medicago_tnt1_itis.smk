# -*- coding: utf-8 -*-
import os
import sys
import time
import pandas as pd
from os.path import join, abspath
from scripts.check_ngspipedb_update import parse_with_wget, parse_with_urllib, ask_if_update_auto
from utils.message import message_success,message_error

# check version
#version = '0.0.3'
#message = parse_with_urllib("http://www.liu-lab.com/ngspipedb/changelog.md", version)
#sys.stderr.write(message)
#print(version)
#ask_if_update_auto()
#time.sleep(5)


shell.executable("bash")

# relative path
snake_dir = workflow.basedir # all configfile, scripts, restructuretext, ens are relative to snakefile (this file)
working_dir = os.getcwd() # input and output path are relative to current working directory

# sub directory
config["resultsDir"] = join(working_dir, config["results_name"], 'ngspipe_result')
config["reportsDir"] = join(working_dir, config["results_name"], 'report')
config["samplesDir"] = config["rawreads_dir"]
config["genomeFasta"] = config["genomeFasta_path"]
config["genomeAnno"] = config["genomeAnno_path"]

# read sample.xls #
smpList = pd.read_csv(config["sample_path"], index_col=0, header=None)
#SAMPLES = list(smpList.index)[0:] # this can be change to [0:2] if you want use the first two sample to run analysis.
SAMPLES = list(smpList.index)[0:] if config['samples_num'] == 'all' else list(smpList.index)[0:int(config['samples_num'])]

# mark run
os.makedirs(join(working_dir, '.ngspipedb'), exist_ok=True)

flag_outdir = join(config["resultsDir"], "flag")
os.makedirs(flag_outdir, exist_ok=True)

# new condition file



# 6. statistic


#report: "report/workflow.rst"
report_outdir = join(config["reportsDir"], "report.html")
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

# final output 1-7
rule all:
    input:
        #
        # 1. sampling data #
        #sampling_data_result                = expand(join(sampling_data_outdir, "{sample}"+config["read1Suffix"]), sample=SAMPLES),
        sampling_reads                      = join(flag_outdir, 'sampling_reads.ok'),
        #
        # 2. raw reads qc #
        rawreads_qc                        = join(flag_outdir, 'rawreads_qc.ok'),
        #
        # 3. tnt1 #
        tnt1                               = join(flag_outdir, 'tnt_merge.ok'),

onsuccess:
    print(message_success)
    shell("python {}/scripts/sendmail.py -r {} -t {} -l {}".format(snake_dir, config['email_addr'], "success", "{log}"))

onerror:
    print(message_error)
    shell("python {}/scripts/sendmail.py -r {} -t {} -l {}".format(snake_dir, config['email_addr'], "error", "{log}"))

include: join("rules", "sampling_reads.smk")
include: join("rules", "rawreads_qc.smk")
include: join("rules", "medicago_tnt1.smk")