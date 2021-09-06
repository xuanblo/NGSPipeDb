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
  'outpdf'       , 'o', 1, "character",'sig_heatmap_Plot outputfile'
), byrow=TRUE, ncol=5)

args = getopt(command)

if ( !is.null(args$help) || is.null(args$deseq2RData) || is.null(args$outpdf) ) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}

library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(DESeq2)

load(args$deseq2RData)
#load("results/rnaseq/diff/diff_by_deseq2/deseq2.RData")

vst <- varianceStabilizingTransformation(dds, blind=FALSE)

# 基因表达量方差过滤
topVarGenes <- head(order(rowVars(assay(vst)), decreasing = TRUE), dim(vst)[1]*0.1) # 取百分之10

### Extract normalized expression for significant genes from the OE and control samples (2:4 and 7:9)
norm_OEsig <- assay(vst[topVarGenes,])

### Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

### Run pheatmap using the metadata data frame for the annotation
pdf(NULL)
pheatmap(norm_OEsig, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = T,
         cluster_cols = T, 
         fontsize_row = 10, 
         border_color = NA, 
         filename = args$outpdf
)
