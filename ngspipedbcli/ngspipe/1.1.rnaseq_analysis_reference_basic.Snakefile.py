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

def list_basename(directory_path, remove_pattern):
    '''
    get file pattern in report directory
    '''
    import glob
    files = glob.glob(directory_path + '/*.pdf')
    basenames = [os.path.basename(i) for i in files]
    return [i.replace(remove_pattern, '') for i in basenames]

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
config["conditionPath"] = config["condition_path"]

# read sample.xls #
smpList = pd.read_csv(config["sample_path"], index_col=0, header=None)
#SAMPLES = list(smpList.index)[0:] # this can be change to [0:2] if you want use the first two sample to run analysis.
SAMPLES = list(smpList.index)[0:] if config['samples_num'] == 'all' else list(smpList.index)[0:int(config['samples_num'])]

# mark run
os.makedirs(join(working_dir, '.ngspipedb'), exist_ok=True)

flag_outdir = join(config["resultsDir"], "flag")
os.makedirs(flag_outdir, exist_ok=True)

# new condition file
condition_df = pd.read_csv(config["condition_path"], index_col=0, header=0)
condition_df = condition_df.loc[SAMPLES, :]
new_conditionpath = join(flag_outdir, 'condition.csv')
condition_df.to_csv(new_conditionpath)


# 6. statistic
statistic_data_all = [
                  '0.genomeFa', 
                  '0.genomeAnno', 
                  '1.rawReads', 
                  '2.cleanReads', 
                  '2.multiqc', 
                  '3.bam', 
                  '4.mergedGtf', 
                  '5.exp',
                  ]
genomeFa_outdir, genomeAnno_outdir, rawReads_outdir, cleanReads_outdir, multiqc_outdir, bam_outdir, mergedGtf_outdir, exp_outdir \
  = [join(config["resultsDir"], "statistic", "statistic_data_of_{}".format(i)) for i in statistic_data_all]

statistic_data_choose = [
                  #'0.genomeFa', 
                  #'0.genomeAnno', 
                  '1.rawReads', 
                  '2.cleanReads', 
                  #'2.multiqc', 
                  '3.bam', 
                  #'4.mergedGtf', 
                  #'5.exp',
                  ]
stat_outdir = join(config["resultsDir"], "statistic")


# 10. network
coexp_outdir = join(config["resultsDir"], "coexp_outdir")

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
        # 3. juntion alignment #
        mapping                            = join(flag_outdir, 'mapping.ok'),
        #
        # 4. assembly
        transcript_assembly                = join(flag_outdir, 'transcript_assembly.ok'),
        #
        # 5. qunatification #
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
        #network = join(coexp_outdir, 'co_exp.network.tsv'),

onstart:
    shell('echo 0.0.19')
    
onsuccess:
    print(message_success)
    shell("python {}/scripts/sendmail.py -r {} -t {} -l {}".format(snake_dir, config['email_addr'], "success", "{log}"))

onerror:
    print(message_error)
    shell("python {}/scripts/sendmail.py -r {} -t {} -l {}".format(snake_dir, config['email_addr'], "error", "{log}"))

include: join("rules", "sampling_reads.smk")
#include: join("rules", "3.junction_align_by_{}.Snakefile.py".format(junction_align_method))
#include: join("rules", "5.quant_by_{}_{}.Snakefile.py".format(quantify_method, 'basic'))
#include: join("rules", "6.statistic_data_of_bam.Snakefile.py")
#include: join("rules", "6.statistic_data_of_rawReads.Snakefile.py")
#include: join("rules", "6.statistic_data_of_cleanReads.Snakefile.py")
#include: join("rules", "7.report_rnaseq_basic.Snakefile.py")
#include: join("rules", "9.differential_expression.deseq2.Snakefile.py")
#include: join("rules", "annotation_protein_by_{}.Snakefile.py".format(anno_method))
#include: join("rules", "differential_expressed_gene_enrich.Snakefile.py")
#include: join("rules", "network_analysis_by_gcen.Snakefile.py")
include: join("rules", "rawreads_qc.smk")
include: join("rules", "transcriptome_mapping.smk")
include: join("rules", "transcript_assembly.smk")
include: join("rules", "expression_quantification.smk")
include: join("rules", "differential_expression.smk")
include: join("rules", "protein_annotation.smk")
include: join("rules", "enrich.smk")