import os

BASE_DIR = os.path.dirname(os.path.dirname(__file__))

GFF_PATH = os.path.join(BASE_DIR, 'PVDB.genome.gff_anno.sqlite3')

SITE_TITLE = 'Plukenetia volubilis'

SERVER_NAME = 'http://pvdb.xtbg.ac.cn'

FTP_NAME = 'ftp://pvdb.xtbg.ac.cn'

GENOME_FASTA = os.path.join(BASE_DIR, 'PVDB.genome.fa')