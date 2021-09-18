import os
import gffutils
from optparse import OptionParser

usage = '''
    python filter_exon_step1.py -i merged.gtf -n exon_number -o step1_exon_filtered.gtf
'''

parser = OptionParser(usage)
parser.add_option('-i', '--input-file', dest='gtfin', help='merged gtf (stringtie merge output)', action='store', type='string')
parser.add_option('-n', '--exon_number', dest='exonNum', help='exon bellow this number will be filter', action='store', type='int', default='2')
parser.add_option('-o', '--output-file', dest='gtfout', help='gtf file after filter exon', action='store', type='string')

(options, args) = parser.parse_args()

os.system('echo Reading gtf file ... && date')

db = gffutils.create_db(options.gtfin, ':memory:')

# transcript / mRNA / lncRNA
transcript_type = 'transcript'

os.system('echo Deleting exon ... && date')

for transcript in db.features_of_type(transcript_type):
    exonCounts = len(list(db.children(transcript, featuretype='exon')))
    if exonCounts < options.exonNum:
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
