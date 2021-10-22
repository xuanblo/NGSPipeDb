#!/usr/env/bin/Rscript
# _*_ coding: utf-8 _*_

#options(getopt.quiet = TRUE)
if (!require('getopt', character.only = TRUE)) {
  install.packages('getopt', dependencies = TRUE, quiet = TRUE)
  #suppressMessages(library('getopt', character.only = TRUE))
  library('getopt', character.only = TRUE)
}

command = matrix(c(
  'help'           , 'h', 0, "logical", 'help message',
  'stringtie_quant_outdir'   , 'i', 1, "character", 'restringtie quant output directory',
  'conditionfile'  , 'c', 1, "character",'flaged sample info',
  'tx2gene'  , 't', 1, "character",'transcript id to gene id',
  'output', 'o', 1, "character",'output file'
), byrow=TRUE, ncol=5)

args = getopt(command)

if ( !is.null(args$help) || is.null(args$stringtie_quant_outdir) || is.null(args$conditionfile) || is.null(args$tx2gene) || is.null(args$output)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}

library(tximport)

samps <- read.table(file.path(args$conditionfile), sep=',', head=TRUE) 
#head(samps)
samps$Sample <- factor(samps$Sample) 
samps$Tissue <- factor(samps$Tissue) 
#table(samps$Sample)
files <- file.path(args$stringtie_quant_outdir, samps$sample_id, "t_data.ctab") 
names(files) <- samps$sample_id 
#head(files)
library(readr)
tmp <- read_tsv(files[1])
tx2gene <- tmp[, c("t_name", "gene_id")]

#head(tx2gene)

txi <- tximport(files, type = "stringtie", tx2gene = tx2gene)

write.table(txi$counts, args$output, sep=',', quote=FALSE)