# add gene/exon feature
agat_sp_gxf_to_gff3.pl -g xingyouteng.gene.final.changeid.gff -o Pvo.gff

python gffrename.py Pvo.gff mRNA_old2new.csv >PVDB.genome.gff

python fastarename.py Function_annotation.xls mRNA_old2new.csv Function_annotation.csv xingyouteng.gene.final.changeid.pep >PVDB.pep.fa
python fastarename.py Function_annotation.xls mRNA_old2new.csv Function_annotation.csv xingyouteng.gene.final.changeid.cds >PVDB.cds.fa

python gffanno.py Function_annotation.xls mRNA_old2new.csv PVDB.genome.gff >PVDB.genome.note.gff

faSize -detailed PVDB.genome.fa >PVDB.genome.size
gff3ToGenePred PVDB.genome.gff PVDB.genome.pred
genePredToBed PVDB.genome.pred PVDB.genome.bed
bedSort PVDB.genome.bed PVDB.genome.sorted.bed
python bedanno.py Function_annotation.xls mRNA_old2new.csv PVDB.genome.sorted.bed >PVDB.genome.sorted.bed15
bedToBigBed -type=bed12+3 -as=bigGenePred.as PVDB.genome.sorted.bed15 PVDB.genome.size PVDB.genome.sorted.bb
faToTwoBit PVDB.genome.fa PVDB.genome.2bit

gffutils-cli create --output PVDB.genome.gff_anno.sqlite3 PVDB.genome.gff
