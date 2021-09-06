#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 10:51:26 2019

@author: zhangxuan
"""

import os
import re
import sys
import pandas as pd

def percentage(num, total):
    p = round(num/total*100,2)
    return str(num) + "\n(" + str(p) + "%)"

files = sys.argv[1:-1]
resultfile = sys.argv[-1]
samples = []
total_reads = []
total_mapped = []
multiple_mapped = []
uniquely_mapped = []
read_1 = []
read_2 = []
read2f = []
read2r = []
nonsplice_reads = []
splice_reads = []
proper_pairs = []
proper_pairs2chr = []


for rseqc in files:
    samplefile = os.path.basename(rseqc)
    sample = samplefile.split(".")[0]
    samples.append(sample)
    print("PROCESS: {}".format(sample))
    with open(rseqc, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                matchObj = re.match(r'(.*):\s*(\d+)', line)
                if matchObj:
                    if matchObj.group(1) == "Total records":
                        total_reads.append(int(matchObj.group(2)))
                    if matchObj.group(1) == "mapq < mapq_cut (non-unique)":
                        multiple_mapped.append(int(matchObj.group(2)))
                    if matchObj.group(1) == "mapq >= mapq_cut (unique)":
                        uniquely_mapped.append(int(matchObj.group(2)))
                    if matchObj.group(1) == "Read-1":
                        read_1.append(int(matchObj.group(2)))
                    if matchObj.group(1) == "Read-2":
                        read_2.append(int(matchObj.group(2)))
                    if matchObj.group(1) == "Reads map to '+'":
                        read2f.append(int(matchObj.group(2)))
                    if matchObj.group(1) == "Reads map to '-'":
                        read2r.append(int(matchObj.group(2)))
                    if matchObj.group(1) == "Non-splice reads":
                        nonsplice_reads.append(int(matchObj.group(2)))
                    if matchObj.group(1) == "Splice reads":
                        splice_reads.append(int(matchObj.group(2)))
                    if matchObj.group(1) == "Reads mapped in proper pairs":
                        proper_pairs.append(int(matchObj.group(2)))
                    if matchObj.group(1) == "Proper-paired reads map to different chrom":
                        proper_pairs2chr.append(int(matchObj.group(2)))
    total_mapped.append(multiple_mapped[-1] + uniquely_mapped[-1])
    
    total_mapped[-1] = percentage(total_mapped[-1], total_reads[-1])
    multiple_mapped[-1] = percentage(multiple_mapped[-1], total_reads[-1])
    uniquely_mapped[-1] = percentage(uniquely_mapped[-1], total_reads[-1])
    read_1[-1] = percentage(read_1[-1], total_reads[-1])
    read_2[-1] = percentage(read_2[-1], total_reads[-1])
    read2f[-1] = percentage(read2f[-1], total_reads[-1])
    read2r[-1] = percentage(read2r[-1], total_reads[-1])
    nonsplice_reads[-1] = percentage(nonsplice_reads[-1], total_reads[-1])
    splice_reads[-1] = percentage(splice_reads[-1], total_reads[-1])
    proper_pairs[-1] = percentage(proper_pairs[-1], total_reads[-1])
    proper_pairs2chr[-1] = percentage(proper_pairs2chr[-1], total_reads[-1])
    
data = {"Sample name": samples,
        "Total reads": total_reads,
        "Total mapped": total_mapped,
        "Multiple mapped": multiple_mapped,
        "Uniquely mapped": uniquely_mapped,
        "Read-1": read_1,
        "Read-2": read_2,
        "Reads map to '+'": read2f,
        "Reads map to '-'": read2r,
        "Non-splice reads": nonsplice_reads,
        "splice reads": splice_reads,
        "Reads mapped in proper pairs": proper_pairs,
        "Proper-paired reads map to different chrom": proper_pairs2chr,
        }

df = pd.DataFrame(data)

df.to_csv(resultfile)