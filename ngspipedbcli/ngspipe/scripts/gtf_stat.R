# gtf statistic

# https://cran.r-project.org/web/packages/refGenome/vignettes/refGenome.pdf
# https://davetang.org/muse/2017/08/04/read-gtf-file-r/

install.packages("refGenome")
library(refGenome)

# create ensemblGenome object for storing Ensembl genomic annotation data
ens <- ensemblGenome()

# read GTF file into ensemblGenome object
read.gtf(ens, "transcript.gtf")

# counts all annotations on each seqname
head(tableSeqids(ens))

# returns all annotations on mitochondria
extractSeqids(ens, 'NW_012126474.1')

# summarise features in GTF file
tableFeatures(ens)

# create table of genes
my_gene <- getGenePositions(ens)
dim(my_gene)

# gene IDs are unique
length(my_gene$gene_id)

library(dplyr)
# use dplyr to create more summaries
# number of genes on each seqname
my_gene %>% group_by(seqid) %>% summarise(n())
my_gene %>% filter(seqid == "NW_012126474.1") %>% select(gene_id) %>% head()
my_gene %>% filter(seqid == 2) %>% select(gene_id) %>% head()

# length of genes
my_gene_length <- my_gene$end - my_gene$start
my_density <- density(my_gene_length)
plot(my_density, main = 'Distribution of gene lengths')

table(my_gene$gene_biotype)