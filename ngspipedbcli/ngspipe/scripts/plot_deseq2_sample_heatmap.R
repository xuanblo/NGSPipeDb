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
  'outpdf'       , 'o', 1, "character",'sig_heatmapPlot outputfile'
), byrow=TRUE, ncol=5)

args = getopt(command)

if ( !is.null(args$help) || is.null(args$deseq2RData) || is.null(args$outpdf) ) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}

library(DESeq2)
library(ggplot2)

load(args$deseq2RData)

vst <- varianceStabilizingTransformation(dds, blind=FALSE)


### Extract the rlog matrix from the object
vst_mat <- assay(vst)    ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2

### Compute pairwise correlation values
vst_cor <- cor(vst_mat)    ## cor() is a base R function

#head(vst_cor)   ## check the output of cor(), make note of the rownames and colnames

### Load pheatmap package
library(pheatmap)

### Plot heatmap
pdf(NULL)
pheatmap(vst_cor, border_color=NA, fontsize = 10, fontsize_row = 10)
pheatmap(vst_cor, border_color=NA, fontsize = 10, fontsize_row = 10, filename=args$outpdf)

