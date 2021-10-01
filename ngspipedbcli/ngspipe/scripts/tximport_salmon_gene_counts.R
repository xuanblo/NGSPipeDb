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
  'salmon_quant_outdir'   , 'i', 1, "character", 'salmon quant output directory',
  'conditionfile'  , 's', 1, "character",'flaged sample info',
  'tx2gene'  , 'g', 1, "character",'tx2gene file',
  'output', 'o', 1, "character",'output file'
), byrow=TRUE, ncol=5)
args = getopt(command)

if ( !is.null(args$help) || is.null(args$salmon_quant_outdir) || is.null(args$conditionfile) || is.null(args$tx2gene) || is.null(args$output)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}

library(tximport)

samps <- read.table(file.path(args$conditionfile), sep=',', head=TRUE) 
#head(samps)
samps$Sample <- factor(samps$Sample) 
samps$Tissue <- factor(samps$Tissue) 
table(samps$Sample)
files <- file.path(args$salmon_quant_outdir, samps$sample_id, "quant.sf") 
names(files) <- samps$sample_id 
#head(files)

tx2gene <- read.table(args$tx2gene, header=FALSE)

#head(tx2gene)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

write.table(txi$counts, args$output, sep=',', quote=FALSE)