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
