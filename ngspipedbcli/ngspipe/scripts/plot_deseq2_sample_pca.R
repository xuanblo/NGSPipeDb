#!/usr/env/bin/Rscript
# _*_ coding: utf-8 _*_

#options(getopt.quiet = TRUE)
if (!require('getopt', character.only = TRUE)) {
  install.packages('getopt', dependencies = TRUE, quiet = TRUE)
  #suppressMessages(library('getopt', character.only = TRUE))
  library('getopt', character.only = TRUE)
}
command = matrix(c(
  'help'          , 'h', 0, "logical", 'help message',
  'deseq2RData'   , 'i', 1, "character", 'saved deseq2 RData',
  'outpdf'       , 'o', 1, "character",'pcaPlot outputfile'
), byrow=TRUE, ncol=5)

args = getopt(command)

if ( !is.null(args$help) || is.null(args$deseq2RData) || is.null(args$outpdf) ) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}

## 设置默认值
#if ( is.null(args$pdfinput) ) {
#  args$pdfinput = 0
#  }

library(DESeq2)
library(ggplot2)


load(args$deseq2RData)
#load("results/rnaseq/diff/diff_by_deseq2/deseq2.RData")

vst <- varianceStabilizingTransformation(dds, blind=FALSE)
 

#BiocManager::install('PCAtools')
library(PCAtools)
pdf(NULL)
p <- pca(assay(vst), metadata = coldata, removeVar = 0.1)

biplot(p,
       colby = colnames(coldata)[1],
       lab=p$yvars,
       hline = 0, vline = c(-25, 0, 25),
       legendPosition = 'top', legendLabSize = 10, legendIconSize = 8.0)

ggsave(filename=args$outpdf, device = "pdf")

