import sys
import os
import datetime
from optparse import OptionParser
from scipy import stats
import matplotlib.pyplot as plt
import gffutils

inputFile = sys.argv[1]
usage = '''
    After assembly, we can get GFF or GTF files. Statistics of 
    annotated files can be used as a reference for adjusting assembly
    parameters.

    The default output folder is gtf2stat.

    The input file can be a GTF files.

    Alternatively, you can specify the filenames directly with
    :option:`-i`/:option:`--input-file` and
    :option:`-o`/:option:`--output-dir` option.

    Example::

        gff_stat_vistualizition.py -i test.gtf -o output_folder
'''

parser = OptionParser(usage)
parser.add_option('-i', '--input-file', dest='gtf', help='read gtf or gff file', metavar="FILE", action = "store", type="string")
parser.add_option("-o","--output-dir", dest="gtf2statdir", help="save file to this directory", default="gtf2statdir", metavar="FILE", action = "store", type="string")

(options, args)=parser.parse_args()

# 读gtf文件
db = gffutils.create_db(options.gtf, ':memory:', merge_strategy='create_unique')

# outputdir
outdir = options.gtf2statdir

if not os.path.isdir(outdir):
    os.makedirs(outdir, exist_ok=True)

os.system('date')

step1 = '''
# ---------------------
# 1. 统计gff中的特征类型
# ---------------------
'''
print(step1)

print("Feature types = ", list(db.featuretypes()))

step2 = '''
# ---------------------
# 2. 统计各个feature数量
# ---------------------
'''
print(step2)

geneSubFeature = "transcript"

geneCounts = db.count_features_of_type(featuretype="gene")
mRNACounts = db.count_features_of_type(featuretype=geneSubFeature)

print({"gene": geneCounts, "transcript": mRNACounts})

## draw plot

plt.figure()
plt.bar(['gene', 'transcript'], [geneCounts, mRNACounts])
plt.title("Feature number barplot")

plt.savefig(os.path.join(outdir, "Feature_number_barplot.pdf"), format='pdf', bbox_inches='tight')

step3 = '''
# --------------------------------------------
# 3. 统计一个基因有几个转录本，一个转录本有多少个exon
# --------------------------------------------
'''
print(step3)

## gene
genes2transcriptCountList = []

for gene in db.features_of_type("gene"):
    transcriptCounts = len(list(db.children(gene, featuretype=geneSubFeature)))
    genes2transcriptCountList.append(transcriptCounts)

print("genes2transcriptCountList:", stats.describe(genes2transcriptCountList))

## exon
transcripts2exonCountList = []

for transcript in db.features_of_type(geneSubFeature):
    exonCounts = len(list(db.children(transcript, featuretype='exon')))
    transcripts2exonCountList.append(exonCounts)

print("transcripts2exonCountList:", stats.describe(transcripts2exonCountList))

## intron
db.update(db.create_introns())

transcripts2intronCountList = []

for transcript in db.features_of_type(geneSubFeature):
    intronCounts = len(list(db.children(transcript, featuretype='intron')))
    if intronCounts:
        transcripts2intronCountList.append(intronCounts)

print("transcripts2intronCountList:", stats.describe(transcripts2intronCountList))

## plot

plt.figure(figsize=(8,8))
plt.suptitle("Feature children number distribution")
plt.subplots_adjust(hspace =0.4)

plt.subplot(3, 1, 1)
plt.hist(genes2transcriptCountList)
plt.title("transcript number per gene")

plt.subplot(3, 1, 2)
plt.hist(transcripts2exonCountList)
plt.title("exon number per transcript")

plt.subplot(3, 1, 3)
plt.hist(transcripts2intronCountList)
plt.title("intron number per transcript")

plt.savefig(os.path.join(outdir, "Feature_children_number_distribution.pdf"), format='pdf', bbox_inches='tight')

step4 = '''
# --------------------------------------------
# 4. 统计各个feature长度分布
# --------------------------------------------
'''
print(step4)

# gene
geneLens = []

for gene in db.features_of_type("gene"):
    gLen = gene.end - gene.start + 1
    geneLens.append(gLen)

print("geneLens:", stats.describe(geneLens))

# transcript
transcriptLens = []

for transcript in db.features_of_type(geneSubFeature):
    tLen = transcript.end - transcript.start + 1
    transcriptLens.append(tLen)

print("transcriptLens:", stats.describe(transcriptLens))

# exon
exonLens = []

for exon in db.features_of_type("exon"):
    eLen = exon.end - exon.start + 1
    exonLens.append(eLen)

print("exonLens:", stats.describe(exonLens))

# intron
intronLens = []

for intron in db.features_of_type("intron"):
    iLen = intron.end - intron.start + 1
    intronLens.append(iLen)

print("intronLens:", stats.describe(intronLens))


# 基因间区长度分布

chroms = [i['seqid'] for i in db.execute('SELECT DISTINCT seqid FROM features;')]

def genes_on_chrom(chroms):
    """
    Yield genes on `chrom`, sorted by start position
    """
    for g in db.features_of_type('gene', order_by='start', limit=(chroms, 0, 1e9)):
        g.strand = '.'
        yield g

def intergenic(chroms):
    """
    Yield intergenic features
    """
    for chrom in chroms:
        genes = genes_on_chrom(chrom)
        for intergenic in db.interfeatures(genes):
            yield intergenic

intergenicLens = []            
            
for intergenic in intergenic(chroms):
    iLen = abs(intergenic.end - intergenic.start) + 1
    intergenicLens.append(iLen)

print("intergenicLens:", stats.describe(intergenicLens))

## plot

plt.figure(figsize=(15,15))
plt.suptitle('Feature length distrubution')
plt.subplots_adjust(hspace =0.4)

plt.subplot(5, 2, 1)
plt.hist(geneLens)
plt.title("gene length")
plt.subplot(5, 2, 2)
plt.boxplot(geneLens)
plt.title("gene length")

plt.subplot(5, 2, 3)
plt.hist(transcriptLens)
plt.title("transcript length")
plt.subplot(5, 2, 4)
plt.boxplot(transcriptLens)
plt.title("transcript length")

plt.subplot(5, 2, 5)
plt.hist(exonLens)
plt.title("exon length")
plt.subplot(5, 2, 6)
plt.boxplot(exonLens)
plt.title("exon length")

plt.subplot(5, 2, 7)
plt.hist(intronLens)
plt.title("intron length")
plt.subplot(5, 2, 8)
plt.boxplot(intronLens)
plt.title("intron length")

plt.subplot(5, 2, 9)
plt.hist(intergenicLens)
plt.title("intergenic length")
plt.subplot(5, 2, 10)
plt.boxplot(intergenicLens)
plt.title("intergenic length")

plt.savefig(os.path.join(outdir, "Feature_length_distrubution.pdf"), format='pdf', bbox_inches='tight')

print()
os.system('date')
print("Finished stat. Three pdfs have been saved.")
