#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 15:44:20 2019

@author: zhangxuan
"""

import os
import sys
import re
import pandas as pd

files = sys.argv[1:-1]
resultfile = sys.argv[-1]

samples = []
clean_reads = []
clean_base = []
error_rate = []
q20 = []
q30 = []
gc = []
error_rate = []
regex = re.compile(r'\s+')

sample_num = int(len(files)/4)

for i in range(int(len(files)/4)):
    read1_basic = files[i+sample_num*0]
    read2_basic = files[i+sample_num*1]
    read1_gc = files[i+sample_num*2]
    read2_gc = files[i+sample_num*3]
    print(i, read1_basic, read2_basic, read1_gc, read2_gc, end='')
    sample_filename1 = os.path.basename(read1_basic)
    sample1 = sample_filename1.split(".")[0]
    sample_filename2 = os.path.basename(read2_basic)
    sample2 = sample_filename2.split(".")[0]
    samples.append(sample1)
    samples.append(sample2)
    print("\tPROCESS: {read1},{read2}".format(read1=sample1,read2=sample2))

    #{sample}_1.seqkit.summary
    with open(read1_basic, 'r') as f:
        for line in f:
            line = line.strip()
            line = line.replace(",", "")
            if "Q30" in line:
                continue
            else:
                items = regex.split(line)
    clean_reads.append(items[3])
    clean_base.append(str(round(int(items[4])/1000000000,2)) + "G")
    q20.append(items[13])
    q30.append(items[14])

    #{sample}/{sample}_1.seqkit.gc_qual
    with open(read1_gc, 'r') as f:
        line = f.readline()
        line = line.strip()
        items = regex.split(line)
    gc.append(round(float(items[0]),2))
    error_rate.append(round(10**(-float(items[1])/10),2))

    #{sample}_2.seqkit.summary
    with open(read2_basic, 'r') as f:
        for line in f:
            line = line.strip()
            line = line.replace(",", "")
            if "Q30" in line:
                continue
            else:
                items = regex.split(line)
    clean_reads.append(items[3])
    clean_base.append(str(round(int(items[4])/1000000000,2)) + "G")
    q20.append(items[13])
    q30.append(items[14])

    #{sample}/{sample}_2.seqkit.gc_qual
    with open(read2_gc, 'r') as f:
        line = f.readline()
        line = line.strip()
        items = regex.split(line)
    gc.append(round(float(items[0]),2))
    error_rate.append(round(10**(-float(items[1])/10),2))

data = {"Sample name": samples,
        "Clean reads": clean_reads,
        "Clean base": clean_base,
        "Error rate(%)": error_rate,
        "Q20(%)": q20,
        "Q30(%)": q30,
        "GC content(%)": gc}

df = pd.DataFrame(data)

df.to_csv(resultfile)

