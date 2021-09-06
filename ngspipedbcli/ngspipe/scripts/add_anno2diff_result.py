import sys
import re
import glob
import os
import csv
import pandas as pd
from urllib.parse import unquote,quote

gtf_file = sys.argv[1]
diff_dir = sys.argv[2]
outfile_dir = sys.argv[3]


if gtf_file.endswith('gtf'):
    gtf_gene_pat = re.compile('gene_id "([^"]+)')
    gtf_transcript_pat = re.compile('gene_id "([^"]+).*transcript_id "([^"]+)')
    gtf_transcript_pat_rev = re.compile('transcript_id "([^"]+).*gene_id "([^"]+)')
else:
    gtf_gene_pat = re.compile('ID=([^;]+)')
    gtf_transcript_pat = re.compile('ID=([^;]+).*Parent=([^;]+)')
    gtf_transcript_pat_rev = re.compile('Parent=([^;]+).*ID=([^;]+)')

gene_anno = dict()

with open(gtf_file, 'r') as f:
    for line in f.readlines():
        if '\tgene\t' in line:
            line = line.rstrip()
            items = line.split('\t')
            matObj = gtf_gene_pat.match(items[8])
            #print(items[8])
            if matObj:
                gid = matObj.group(1)
            else:
                sys.stderr.write('ERROR format in {}\n'.format(items[8]))
                sys.exit(-1)
            gene_anno[gid] = unquote(items[8])
        if '\ttranscript\t' in line or '\tmiRNA\t' in line or '\tmRNA\t' in line or '\tncRNA\t' in line or '\trRNA\t' in line or '\ttRNA\t' in line:
            line = line.rstrip()
            items = line.split('\t')
            matObj = gtf_transcript_pat.match(items[8])
            matObj_rev = gtf_transcript_pat_rev.match(items[8])
            #print(items[8])
            if matObj:
                if gtf_file.endswith('gtf'):
                    gid = matObj.group(1)
                    tid = matObj.group(2)
                else:
                    gid = matObj.group(2)
                    tid = matObj.group(1)
            elif matObj_rev:
                if gtf_file.endswith('gtf'):
                    gid = matObj.group(2)
                    tid = matObj.group(1)
                else:
                    gid = matObj.group(1)
                    tid = matObj.group(2)
            else:
                sys.stderr.write('ERROR format in {}\n'.format(items[8]))
                sys.exit(-1)
            gene_anno[gid] = gene_anno[gid] + ' ^S^ ' + unquote(items[8])

path = os.path.join(diff_dir, '*.csv')

for diff_res in glob.glob(path):
    filename = os.path.basename(diff_res)
    sys.stderr.write('READ: ' + filename + '\n')
    outfile = os.path.join(outfile_dir, filename)

    gene_anno_df = pd.DataFrame.from_dict(gene_anno, orient='index', columns=['anno'])
    
    diff_df = pd.read_csv(diff_res, header=0, index_col=0, sep=',')
    diff_df.index = [i.split('|')[0] for i in diff_df.index]

    diff_df['gene_anno'] = gene_anno_df.loc[diff_df.index,'anno']
    print(diff_df)
    diff_df.to_csv(outfile,quoting=csv.QUOTE_NONE,sep=',',escapechar='|')