import os
import gffutils
from optparse import OptionParser

usage = '''
    python filter_transcript_length_step2.py -i step1_exon_filter.gtf -n trascript_length -o step2_length_filter.gtf
'''

# 过滤基因最长转录本小于200的基因

parser = OptionParser(usage)
parser.add_option('-i', '--input-file', dest='gtfin', help='merged gtf (stringtie merge output)', action='store', type='string')
parser.add_option('-l', '--trascript_length', dest='len', help='transcript bellow this length will be filter', action='store', type='int', default='150')
parser.add_option('-o', '--output-file', dest='gtfout', help='gtf file after filter length', action='store', type='string')

(options, args) = parser.parse_args()

os.system('echo Reading gtf file ... && date')

db = gffutils.create_db(options.gtfin, ':memory:')

# transcript / mRNA / lncRNA
transcript_type = 'transcript'

os.system('echo Deleting samll length transcript ... && date')

for transcript in db.features_of_type(transcript_type):
    transcript_len = transcript.end - transcript.start + 1
    if transcript_len < options.len:
        exons = db.children(transcript)
        db.delete(exons)

for gene in db.features_of_type('gene'):
    db.delete(gene)

for transcript in db.features_of_type(transcript_type):
    db.delete(transcript)

outputfile = options.gtfout

recs = db.all_features()

gff_out = gffutils.gffwriter.GFFWriter(outputfile, in_place=True, with_header=False)

os.system('echo Wrinting ... && date')

gff_out.write_recs(recs)
gff_out.close()
(base) [zhangxuan@t640 RunAll]$ cat filter_low_expression_transcript_step3.py
cat: filter_low_expression_transcript_step3.py: 没有那个文件或目录
(base) [zhangxuan@t640 RunAll]$ cat script/filter_low_expression_transcript_step3.py
import os
import gffutils
import pandas as pd
from optparse import OptionParser

usage = '''
    python filter_low_expression_transcript_step3.py -i step2_length_filter.gtf -n fpkm -o step3_exp_filter.gtf
'''

parser = OptionParser(usage)
parser.add_option('-i', '--input-gtf', dest='gtfin', help='merged gtf (stringtie merge output)', action='store', type='string')
parser.add_option('-m', '--input-exp-matrix', dest='expin', help='expression matrix (fpkm)', action='store', type='string')
parser.add_option('-f', '--fpkm', dest='minfpkm', help='at least in one sample fpkm should be higher than this value', action='store', type='float', default='0.1')
parser.add_option('-o', '--output-file', dest='gtfout', help='gtf file after filter exon', action='store', type='string')

(options, args) = parser.parse_args()

os.system('echo Reading gtf file ... && date')

db = gffutils.create_db(options.gtfin, ':memory:')

# transcript / mRNA / lncRNA
transcript_type = 'transcript'

os.system('echo Reading expression matrix ... && date')

exp_df = pd.read_csv(options.expin, header=0, index_col=0, sep='\t')

max_exp_series = exp_df.max(axis=1)

filter_exp_series = max_exp_series[max_exp_series < options.minfpkm]

for transcript in filter_exp_series.index:
    exons = db.children(transcript)
    db.delete(exons)

for gene in db.features_of_type('gene'):
    db.delete(gene)

for transcript in db.features_of_type(transcript_type):
    db.delete(transcript)

outputfile = options.gtfout

recs = db.all_features()

gff_out = gffutils.gffwriter.GFFWriter(outputfile, in_place=True, with_header=False)

os.system('echo Wrinting ... && date')

gff_out.write_recs(recs)
gff_out.close()
