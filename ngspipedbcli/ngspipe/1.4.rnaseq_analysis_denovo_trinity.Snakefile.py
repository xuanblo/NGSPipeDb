# -*- coding: utf-8 -*-
import os
import sys
import time
import pandas as pd
from os.path import join, abspath
from scripts.check_ngspipedb_update import parse_with_wget, parse_with_urllib, ask_if_update_auto
from utils.message import message_success,message_error

shell.executable("bash")

# relative path
snake_dir = workflow.basedir # all configfile, scripts, restructuretext, ens are relative to snakefile (this file)
working_dir = os.getcwd() # input and output path are relative to current working directory

# sub directory
config["resultsDir"] = join(working_dir, config["results_name"], 'ngspipe_result')
config["reportsDir"] = join(working_dir, config["results_name"], 'report')
config["samplesDir"] = config["rawreads_dir"]
config["conditionPath"] = config["condition_path"]

# read sample.xls #
smpList = pd.read_csv(config["sample_path"], index_col=0, header=None)
SAMPLES = list(smpList.index)[0:] if config['samples_num'] == 'all' else list(smpList.index)[0:int(config['samples_num'])]

# mark run
os.makedirs(join(working_dir, '.ngspipedb'), exist_ok=True)

flag_outdir = join(config["resultsDir"], "flag")
os.makedirs(flag_outdir, exist_ok=True)

# new condition file

condition_df = pd.read_csv(config["condition_path"], index_col=0, header=0)
condition_df = condition_df.loc[SAMPLES, :]
new_conditionpath = join(flag_outdir, 'condition.csv')
if not os.path.exists(new_conditionpath):
    condition_df.to_csv(new_conditionpath)


rule all:
    input:
        #
        # 1. sampling data #
        sampling_reads                      = join(flag_outdir, 'sampling_reads.ok'),
        #
        # 2. raw reads qc #
        rawreads_qc                        = join(flag_outdir, 'rawreads_qc.ok'),
        #
        # 3. assembly
        #trinity_assembly_result                     = join(transcript_assembly_outdir, 'Trinity.fasta'),
        transcript_fa                       = join(transcript_assembly_outdir, 'trinity_assembly2', 'Trinity.samples.tsv'),
        prepareok                           = join(transcript_assembly_outdir, 'prepare_the_reference_for_alignment_and_abundance_estimation', "prepare.ok"),
        trinity_abundance_result            = expand(join(transcript_assembly_outdir, 'Estimating_transcript_abundance', '{sample}', "{}.isoforms.results".format('RSEM')), sample=SAMPLES),
        supertranscript_fa                  = join(transcript_assembly_outdir, 'SuperTranscripts', 'trinity_genes.fasta'),
        #
        # quantification
        quantification                           = join(flag_outdir, "quantification.ok"),
        #
        # 6. differential expression
        differential_expression                      = join(flag_outdir, "differential_expression.ok"),
        #
        # 7. protein annotation
        protein_annotation                       = join(flag_outdir, 'protein_annotation.ok'),
        #
        # 8. go kegg enrich
        enrich                                   = join(flag_outdir, 'enrich.ok'),
        #
        # 9. network
        network                                 = join(flag_outdir, 'network.ok'),

onsuccess:
    print(message_success)
    shell("python {}/scripts/sendmail.py -r {} -t {} -l {}".format(snake_dir, config['emailaddr'], "success", "{log}"))

onerror:
    print(message_error)
    shell("python {}/scripts/sendmail.py -r {} -t {} -l {}".format(snake_dir, config['emailaddr'], "error", "{log}"))


include: join("rules", "rawreads_qc.smk")
include: join("rules", "trinity.smk")
include: join("rules", "expression_quantification.smk")
include: join("rules", "differential_expression.smk")
include: join("rules", "protein_annotation.smk")
include: join("rules", "enrich.smk")
include: join("rules", "network_analysis.smk")

