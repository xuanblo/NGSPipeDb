# -*- coding: utf-8 -*-
import os
from os.path import join
import sys
import pandas as pd
from utils.message import message_success,message_error

# relative path
snake_dir = workflow.basedir # all configfile, scripts, restructuretext, ens are relative to snakefile (this file)
working_dir = os.getcwd() # input and output path are relative to current working directory

# sub directory
config["resultsDir"] = join(working_dir, config["results_name"], 'ngspipe_result')
config["reportsDir"] = join(working_dir, config["results_name"], 'report')
config["samplesDir"] = config["rawreads_dir"]
config["genomeFasta"] = config["genomeFasta_path"]
config["genomeAnno"] = config["genomeAnno_path"]

# ----------------------------------------------------------------------- #
# sample information #
#
smpList = pd.read_csv(config["sample_path"], index_col=0, header=None)
SAMPLES = list(smpList.index)[0:]
# ----------------------------------------------------------------------- #


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------- #
# detail parameters in pipe #
#
# 1. sampling data
# for test the pipe, you can choose to the part of the input file, can be whole,head:40000,tail:40000,random:0.5,random:40000
sampling_method = 'links' # tail, seqkit_number, seqkit_proportion, head, tail
sampling_data_outdir = join(config["resultsDir"], "sampling_data", "sampling_data_by_{}".format(sampling_method))

# 2. raw reads qc
qc_method = 'trim-galore' # trimomatic
qc_outdir = join(config["resultsDir"], "rawReads_qc", "rawReads_qc_by_{}".format(qc_method))

# 3. mapping to genome
mapping_method = 'bwa' # bowtie, bwa, 
mapping_outdir = join(config["resultsDir"], "mapping", "mapping_by_{}".format(mapping_method))
junction_align_outdir = mapping_outdir
genome_index_prefix = 'genome'

# 4. call variation
call_variation_method = 'bcftools' # samtools, bcftools or gatk
variation_outdir = join(config["resultsDir"], "call_variation", "call_by_{}".format(call_variation_method))

# 5. variation annotation
variation_annotation_method = 'snpeff' # snpeff
annotation_outdir = join(config["resultsDir"], "variation_annotation", "annotation_by_{}".format('variation_annotation_method'))
annoFormat = 'gtf22' # gff3, gff2, gtf22
dbname = 'ngspipe-resequencing',

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
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

# ------------------------------
#report: "report/workflow.rst"
report_outdir = join(config["reportsDir"], "report.html")
# ------------------------------

rule all:
    input:
        #
        # 1. sampling data #
        sampling_data_result                = expand(join(sampling_data_outdir, "{sample}"+config["read1Suffix"]), sample=SAMPLES),
        #
        # 2. raw reads qc #
        rawReads_qc_result                 = expand(join(qc_outdir, "{sample}", "{sample}.cleanR1.fq.gz"), sample=SAMPLES),
        #
        # 3. mapping #
        mapping_result                     = expand(join(mapping_outdir, "{sample}", "{sample}.sorted.bam"), sample=SAMPLES),
        #
        # 4. call variation #
        variation_result                   = expand(join(variation_outdir, "{sample}", "{sample}.vcf"), sample=SAMPLES),
        #
        # 5. annotation #
        #
        annotation_result                  = expand(join(annotation_outdir, "{sample}", "{sample}.snpeff.vcf"), sample=SAMPLES),
        #
        # 5. statistic #
        # if raw data doesn't exists, "Failed to solve scheduling problem with ILP solver" error will ocurrs.
        statistic_result                   = expand(join(stat_outdir, "statistic_data_of_{statistic_data}", 'statistic.completed'), statistic_data=statistic_data_choose),
        # 6. report #
        report_result    = join(config['reportsDir'], "report.ok"),

onsuccess:
    print(message_success)
    shell("python {}/scripts/sendmail.py -r {} -t {} -l {}".format(snake_dir, config['email_addr'], "success", "{log}"))

onerror:
    print(message_error)
    shell("python {}/scripts/sendmail.py -r {} -t {} -l {}".format(snake_dir, config['email_addr'], "error", "{log}"))

include: join("rules", "1.sampling_data_by_{}.Snakefile.py".format(sampling_method))
include: join("rules", "2.rawReads_qc_by_{}.Snakefile.py".format(qc_method))
include: join("rules", "3.genome_align_by_{}.Snakefile.py".format(mapping_method))
include: join("rules", "10.call_variation_by_{}.Snakefile.py".format(call_variation_method))
include: join("rules", "11.variation_annotation_by_{}.Snakefile.py".format(variation_annotation_method))
include: join("rules", "6.statistic_data_of_bam.Snakefile.py")
include: join("rules", "6.statistic_data_of_rawReads.Snakefile.py")
include: join("rules", "6.statistic_data_of_cleanReads.Snakefile.py")
include: join("rules", "7.1.report_resequcing.Snakefile.py")