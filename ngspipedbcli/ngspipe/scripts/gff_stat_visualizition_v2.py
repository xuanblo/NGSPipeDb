import sys
import os
import datetime
from optparse import OptionParser
import gffutils
import numpy as np
import pyfaidx
from Bio.SeqUtils import GC
from Bio.SeqUtils import GC123

usage = '''
    After assembly, we can get GFF or GTF files. Statistics of 
    annotated files can be used as a reference for adjusting assembly
    parameters.

    The default output folder is gtf2stat.
        1. summary
        2. distribution file, you can draw plot yourself.

    The input file can be a GTF files.

    Alternatively, you can specify the filenames directly with
    :option:`-i`/:option:`--input-file` and
    :option:`-o`/:option:`--output-dir` option.

    Example::

        gff_stat_vistualizition.py -i test.gtf -o output_folder
'''

parser = OptionParser(usage)
parser.add_option('-i', '--input-file', dest='gtf', help='read gtf or gff file', metavar="FILE", action = "store", type="string")
parser.add_option('-g', '--genome-file', dest='genome', help='read genome file', metavar="FILE", action = "store", type="string")
parser.add_option("-o","--output-dir", dest="gtf2statdir", help="save file to this directory", default="gtf2statdir", metavar="FILE", action = "store", type="string")

os.system('date')

(options, args)=parser.parse_args()

# outputdir
outdir = options.gtf2statdir

if not os.path.isdir(outdir):
    os.makedirs(outdir, exist_ok=True)

# 读gtf文件
db = gffutils.create_db(options.gtf, dbfn=":memory:", merge_strategy="create_unique")
# 读基因组文件
fasta = pyfaidx.Fasta(options.genome)

# report type
items = ["Number of scaffolds",
         "Number of genes",
         "Number of transcripts",
         "Number of exons",
         "Number of intron in exon",
         "Number of single exon gene",
         "Number of single exon transcript",
         "mean transcripts per gene",
         "mean exons per transcript",
         "mean introns in exons per transcript",
         "Total gene length",
         "Total transcript length",
         "Total exon length",
         "Total intron length per exon",
         "mean gene length",
         "mean transcript length",
         "mean exon length",
         "mean intron in exon length",
         "Longest genes",
         "Longest transcripts",
         "Longest exons",
         "Longest intron into exon part",
         "Shortest genes",
         "Shortest transcripts",
         "Shortest exons",
         "Shortest intron into exon part",
         ]

num_scf, num_gene, num_transcript, num_exon, num_intron, num_sgexon_gene, num_sgexon_transcript, mean_transcript_per_gene, mean_exon_per_transcript, mean_intron_per_transcript, total_gene_length, toal_transcript_length, total_exon_length, total_intron_length, mean_gene_length, mean_transcript_length, mean_exon_length, mean_intron_length, longest_gene_num, longest_transcript_num, longest_exon_num, longest_intron_num, shortest_gene_num, shortest_transcript_num, shortest_exon_num, shortest_intron_num, *zeros = np.zeros(100)

arrGeneId = []
arrGeneGC = []
arrGeneLen = []
arrIsoformNum = []
arrTransId = []
arrTransGC = []
arrTransLen = []
arrTransNum = []
arrExonId = []
arrExonGC = []
arrExonLen = []
arrExonNum = []
arrIntronId = []
arrIntronGC = []
arrIntronLen = []
arrIntergenicId = []
arrIntergenicGC = []
arrIntergenicLen = []

geneDistribution = {
    "gene id": arrGeneId,
    "gc": arrGeneGC,
    "length": arrGeneLen,
    "isoform num": arrIsoformNum
}

transcriptDistribution = {
    "transcript id": arrTransId,
    "gc": arrTransGC,
    "length": arrTransLen,
    "exon num": arrTransNum
}

exonDistribution = {
    "exon id": arrExonId,
    "gc": arrExonGC,
    "length": arrExonLen,
    "exon num": arrExonNum
}

intronDistribution = {
    "intron id": arrIntronId,
    "gc": arrIntronGC,
    "length": arrIntronLen,
}

intergenicDistribution = {
    "intergenic id": arrIntergenicId,
    "gc": arrIntergenicGC,
    "length": arrIntergenicLen,
}



# 取序列
#a = genomeObj[4:12]

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

geneSubFeature = "mRNA"

for feature_type in db.featuretypes():
    counts = db.count_features_of_type(featuretype=feature_type)
    print(feature_type, counts)

step3 = '''
# --------------------------------------------
# 3. 统计基因和基因间区不同信息的分布情况
# --------------------------------------------
'''
print(step3)

levelTwoType = []

geneDist = os.path.join(outdir, "gene.dist.tsv")
with open(geneDist, 'w') as f:
    f.write("id\ttranscript_num\tgc\tgc1\tgc2\tgc3\tlength\n")
    for gene in db.features_of_type("gene"):
        transcriptCounts = str(len(list(db.children(gene))))
        transcriptType = [t.featuretype for t in db.children(gene, level=1)]
        levelTwoType += transcriptType
        geneFa = gene.sequence(fasta)
        gc = GC(geneFa)
        gc123 = GC123(geneFa)
        geneLen = gene.end - gene.start + 1
        items = [gene.id, transcriptCounts, str(gc), str(gc123[1]), str(gc123[2]), str(gc123[3]), str(geneLen)]
        linestr = '\t'.join(items)
        f.write(linestr + '\n')


print(set(levelTwoType))

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
        for intergenic in db.interfeatures(features = genes, new_featuretype="intergenic"):
            yield intergenic


intergenicDist = os.path.join(outdir, "intergenic.dist.tsv")
with open(intergenicDist, 'w') as f:
    f.write("id\tgc\tgc1\tgc2\tgc3\tlength\n")
    for intergenic in intergenic(chroms):
        if intergenic.end - intergenic.start < 10:
            continue
        intergenic_id = '-'.join(intergenic.attributes['ID'])
        sequence = intergenic.sequence(fasta)
        gc = GC(sequence)
        gc123 = GC123(sequence)
        intergenicLen = intergenic.end - intergenic.start + 1
        items = [intergenic_id, str(gc), str(gc123[1]), str(gc123[2]), str(gc123[3]), str(intergenicLen)]
        linestr = '\t'.join(items)
        f.write(linestr + '\n')

for type in list(set(levelTwoType)):
    rnaDist = os.path.join(outdir, "{}.dist.tsv".format(type))
    exonDist = os.path.join(outdir, "{}.exon.dist.tsv".format(type))
    intronDist = os.path.join(outdir, "{}.intron.dist.tsv".format(type))
    cdsDist = os.path.join(outdir, "{}.cds.dist.tsv".format(type))
    rna_f = open(rnaDist, 'w')
    exon_f = open(exonDist, 'w')
    intron_f = open(intronDist, 'w')
    cds_f = open(cdsDist, 'w')
    rna_f.write("id\texon_num\tgc\tgc1\tgc2\tgc3\tlength\n")
    exon_f.write("id\tgc\tgc1\tgc2\tgc3\tlength\n")
    intron_f.write("id\tgc\tgc1\tgc2\tgc3\tlength\n")
    cds_f.write("id\tgc\tgc1\tgc2\tgc3\tlength\n")
    for gene in db.features_of_type('gene'):
        for transcript in db.children(gene, featuretype=type):
            exonCounts = str(len(list(db.children(transcript))))
            sequence = transcript.sequence(fasta)
            if not sequence:
                    sys.stderr.write(transcript.id + 'failed\n')
                    continue
            gc = GC(sequence)
            gc123 = GC123(sequence)
            transcriptLen = transcript.end - transcript.start + 1
            items = [transcript.id, str(gc), str(gc123[1]), str(gc123[2]), str(gc123[3]), str(intergenicLen)]
            linestr = '\t'.join(items)
            rna_f.write(linestr + '\n')
            for exon in db.children(transcript, featuretype='exon'):
                sequence = exon.sequence(fasta)
                if not sequence:
                    sys.stderr.write(exon.id + 'failed\n')
                    continue
                gc = GC(sequence)
                gc123 = GC123(sequence)
                exonLen = exon.end - exon.start + 1
                items = [exon.id, exonCounts, str(gc), str(gc123[1]), str(gc123[2]), str(gc123[3]), str(exonLen)]
                linestr = '\t'.join(items)
                exon_f.write(linestr + '\n')
            for intron in db.children(transcript, featuretype='intron'):
                sequence = intron.sequence(fasta)
                if not sequence:
                    sys.stderr.write(intron.id + 'failed\n')
                    continue
                gc = GC(sequence)
                gc123 = GC123(sequence)
                intronLen = intron.end - intron.start + 1
                items = [intron.id, str(gc), str(gc123[1]), str(gc123[2]), str(gc123[3]), str(intronLen)]
                linestr = '\t'.join(items)
                intron_f.write(linestr + '\n')
            for cds in db.children(transcript, featuretype='CDS'):
                sequence = cds.sequence(fasta)
                if not sequence:
                    sys.stderr.write(cds.id + 'failed\n')
                    continue
                gc = GC(sequence)
                gc123 = GC123(sequence)
                cdsLen = cds.end - cds.start + 1
                items = [cds.id, str(gc), str(gc123[1]), str(gc123[2]), str(gc123[3]), str(cdsLen)]
                linestr = '\t'.join(items)
                cds_f.write(linestr + '\n')
    
    rna_f.close()
    exon_f.close()
    intron_f.close()
    cds_f.close()