#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 10:38:39 2019

@author: zhangxuan

Current version 1.0

"""

import os
import re
import sys
import datetime
import pandas as pd
from optparse import OptionParser

usage = '''
    Expression metrics can by FPKM and TPM.

    The result dirs contains stringtie output directory.

    The transcript_matrix_file.

	The gene_matrix_file.

    Alternatively, you can specify the filenames directly with
    :option:`-m`/:option:`--expression_metric_type` and
    :option:`-s`/:option:`--stringtie_result_dirs` 
	:option:`-t`/:option:`--transcript_matrix_file` option.
	:option:`-g`/:option:`--gene_matrix_file` option.

    Example::

        gff_stat_vistualizition.py -m expression_metric -s result_dirs -t transcript_matrix_file -g gene_matrix_file
'''

parser = OptionParser(usage)
parser.add_option('-m', '--expression_metric_type', dest='type', help='fpkm or tpm', metavar="STR", action = "store", type="string")
parser.add_option("-s","--stringtie_result_dirs", dest="stringtiedir", help="a directory with all sample's stringtie output file", metavar="DIR", action = "store", type="string")
parser.add_option("-t","--transcript_matrix_file", dest="tmatrix", help="transcripts expression matrix", metavar="FILE", action = "store", type="string")
parser.add_option("-g","--gene_matrix_file", dest="gmatrix", help="genes expression matrix", metavar="FILE", action = "store", type="string")

(options, args)=parser.parse_args()



def parse_transcript(basedir, sample, type):
    global transcripts
    global tdict
    transcript_file = os.path.join(basedir, sample, sample+".gtf")
    with open(transcript_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.rstrip()
            items = line.split('\t')
            if items[2] == "transcript":
                matchObj = re.match(r'gene_id "(\S+)"; transcript_id "(\S+)";.*cov "(\S+)"; FPKM "(\S+)"; TPM "(\S+)";', items[8])
                if matchObj:
                    gene_id = matchObj.group(1)
                    transcript_id = matchObj.group(2)
                    fpkm = matchObj.group(4)
                    tpm = matchObj.group(5)
                    tdict.setdefault(transcript_id, {})
                    transcripts[transcript_id] = ''
                    if type == 'fpkm':
                        tdict[transcript_id][sample] = fpkm
                    elif type == 'tpm':
                        tdict[transcript_id][sample] = tpm
                    else:
                        sys.exit("not such expr type")
                else:
                    sys.exit("gtf file format not good")

def parse_gene(basedir, sample, type):
    global genes
    global gdict
    gene_file = os.path.join(basedir, sample, sample+'.tab')
    with open(gene_file, 'r') as f:
        for line in f:
            if line.startswith('Gene'):
                continue
            line = line.rstrip()
            items = line.split('\t')
            gene_id = items[0]
            fpkm = items[7]
            tpm = items[8]
            genes[gene_id] = ''
            gdict.setdefault(gene_id, {})
            if type == 'fpkm':
                gdict[gene_id][sample] = fpkm
            elif type == 'tpm':
                gdict[gene_id][sample] = tpm
            else:
                sys.exit("not such expr type")


basedir = options.stringtiedir

samples = [d for d in os.listdir(basedir) if os.path.isdir(os.path.join(basedir, d))]

type = options.type

transcripts = dict()
genes = dict()
tdict = dict()
gdict = dict()


for sample in samples:
    sample_path = os.path.join(basedir, sample)
    if os.path.isfile(sample_path):
        continue
    #Parse gene expression values
    parse_gene(basedir, sample, type)

    #Parse transcript expression values
    parse_transcript(basedir, sample, type)


#Write out the transcript file
tout = open(options.tmatrix, 'w')
tout.write("Transcript_id\t" + '\t'.join(samples) + "\n")
for transcript_id in transcripts:
    line = []
    for sample in samples:
        if tdict[transcript_id][sample]:
            line.append(tdict[transcript_id][sample])
        else:
            line.append(0)
    linestr = [str(i) for i in line]
    tout.write(transcript_id + "\t" + '\t'.join(linestr) + "\n")

#Write out the gene file
gout = open(options.gmatrix, 'w')
gout.write("Gene_id\t" + '\t'.join(samples) + "\n")
for gene_id in genes:
    line = []
    for sample in samples:
        if gdict[gene_id][sample]:
            line.append(gdict[gene_id][sample])
        else:
            line.append(0)
    linestr = [str(i) for i in line]
    gout.write(gene_id + "\t" + '\t'.join(linestr) + "\n")

tout.close()
gout.close()
