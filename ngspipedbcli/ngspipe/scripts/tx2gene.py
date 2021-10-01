import sys
import re
import glob
import os
import csv
import pandas as pd
from urllib.parse import unquote,quote

gtf_file = sys.argv[1]

gtf_gene_pat = re.compile('gene_id "([^"]+)')
gtf_transcript_pat = re.compile('gene_id "([^"]+).*transcript_id "([^"]+)')
gtf_transcript_pat_rev = re.compile('transcript_id "([^"]+).*gene_id "([^"]+)')


tx2gene = dict()

with open(gtf_file, 'r') as f:
    for line in f.readlines():
        line = line.rstrip()
        items = line.split('\t')

        matObj = gtf_transcript_pat.match(items[8])
        matObj_rev = gtf_transcript_pat_rev.match(items[8])

        if matObj:
            gid = matObj.group(1)
            tid = matObj.group(2)
        elif matObj_rev:
            gid = matObj_rev.group(2)
            tid = matObj_rev.group(1)
        else:
            sys.stderr.write('ERROR format in {}\n'.format(items[8]))
            sys.exit(-1)
        tx2gene[tid] = gid

for t,g in tx2gene.items():
    print(t, g, sep='\t')

