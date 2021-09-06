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
#ask_if_update_auto()
#time.sleep(5)


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
SAMPLES = list(smpList.index)[0:] # this can be change to [0:2] if you want use the first two sample to run analysis.

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------- #
# detail parameters in rna-seq analysis #
#
# 1. sampling data
# for test the pipe, you can choose to the part of the input file, can be whole,head:40000,tail:40000,random:0.5,random:40000
sampling_method = 'links' # tail, seqkit_number, seqkit_proportion, head, tail
sampling_data_outdir = join(config["resultsDir"], "sampling_data", "sampling_data_by_{}".format(sampling_method))

# 2. raw reads qc
qc_method = 'trim-galore' # trimmomatic or cutadate or fastp or fastqc or multiqc
qc_outdir = join(config["resultsDir"], "rawReads_qc", "rawReads_qc_by_{}".format(qc_method))

# 3. junction alignmnet
junction_align_method = 'hisat2' # star
junction_align_outdir = join(config["resultsDir"], "mapping", "junction_align_by_{}".format(junction_align_method))
genome_index_prefix = "genome"
# rna-seq sequencing type, can be fr-firststrand, none, fr-secondstrand
rna_library = "" # "--rna-strandness RF"(fr-firststrand) or "--rna-strandness FR"(fr-secondstrand)

# 4. transcript assembly
transcript_assembly_method = 'stringtie' # star
stringtie_rna_library = "" # --rf: Assumes a stranded library fr-firststrand. --fr: Assumes a stranded library fr-secondstrand.
transcript_assembly_outdir = join(config["resultsDir"], "transcript_assembly", "transcript_assembly_by_{}".format(transcript_assembly_method))

# 5. quantification
quantify_method = 'stringtie' # htseqcounts or featurecounts
quantify_outdir = join(config["resultsDir"], "quantify", "quantify_by_{}".format(quantify_method))

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

# 7. diff gene discovery
diff_outdir = join(config["resultsDir"], "diff", "diff_by_{}".format('deseq2'))
diff_anno_outdir = join(config["resultsDir"], "diff", "diff_by_{}".format('deseq2'), "add_gene_note")

#report: "report/workflow.rst"
report_outdir = join(config["reportsDir"], "report.html")
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

# final output 1-7
rule all:
    input:
        #
        # 1. sampling data #
        sampling_data_result                = expand(join(sampling_data_outdir, "{sample}"+config["read1Suffix"]), sample=SAMPLES),
        #
        # 2. raw reads qc #
        rawReads_qc_result                 = expand(join(qc_outdir, "{sample}", "{sample}.cleanR1.fq.gz"), sample=SAMPLES),
        #
        # 3. juntion alignment #
        junction_align_result              = expand(join(junction_align_outdir, "{sample}", "{sample}.sorted.bam"), sample=SAMPLES),
        #
        # 4. qunatification #
        quantify                           = join(quantify_outdir, "gene.csv"),
        #
        # 5. statistic #
        # if raw data doesn't exists, "Failed to solve scheduling problem with ILP solver" error will ocurrs.
        statistic_result                   = expand(join(stat_outdir, "statistic_data_of_{statistic_data}", 'statistic.completed'), statistic_data=statistic_data_choose),
        #
        # 6. report #
        report_result                      = join(config['reportsDir'], "report.ok"),
        #
        # 7. differential expression
        diff_outputok                      = join(diff_anno_outdir, "diff.ok"),
    output:
        touch(join(working_dir, '.ngspipedb', 'run_all.runned'))


onsuccess:
    print(message_success)
    shell("python {}/scripts/sendmail.py -r {} -t {} -l {}".format(snake_dir, config['email_addr'], "success", "{log}"))

onerror:
    print(message_error)
    shell("python {}/scripts/sendmail.py -r {} -t {} -l {}".format(snake_dir, config['email_addr'], "error", "{log}"))

include: join("rules", "1.sampling_data_by_{}.Snakefile.py".format(sampling_method))
include: join("rules", "2.rawReads_qc_by_{}.Snakefile.py".format(qc_method))
include: join("rules", "3.junction_align_by_{}.Snakefile.py".format(junction_align_method))
include: join("rules", "5.quant_by_{}_{}.Snakefile.py".format(quantify_method, 'basic'))
include: join("rules", "6.statistic_data_of_bam.Snakefile.py")
include: join("rules", "6.statistic_data_of_rawReads.Snakefile.py")
include: join("rules", "6.statistic_data_of_cleanReads.Snakefile.py")
include: join("rules", "7.report_rnaseq_basic.Snakefile.py")
include: join("rules", "9.differential_expression.deseq2.Snakefile.py")
