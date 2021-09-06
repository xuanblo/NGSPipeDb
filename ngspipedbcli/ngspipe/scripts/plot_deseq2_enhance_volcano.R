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
  'outdir'       , 'o', 1, "character",'volcano outputfiledir'
), byrow=TRUE, ncol=5)

args = getopt(command)

if ( !is.null(args$help) || is.null(args$deseq2RData) || is.null(args$outdir) ) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}

#BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
library("DESeq2")

load(args$deseq2RData)

diff_result_dir = paste(args$outdir , '/', sep='')

sample_comb <- combn(unique(coldata$Sample),2)

for (i in 1:dim(sample_comb)[2]){
  sample1 <- as.character(sample_comb[, i][1])
  sample2 <- as.character(sample_comb[, i][2])
  
  res <- results(dds, contrast=c("Sample",sample1,sample2))
  #res <- subset(res, padj < 0.01 & abs(log2FoldChange) > 2.32)
  
  # 火山图,在进行MA绘图之前，我们使用 lfcShrink函数缩小log2倍变化
  #res.shrink <- lfcShrink(dds, contrast=c("condition",sample1,sample2), res=res)
  #dev.new()
  dir.create(diff_result_dir, recursive = TRUE)
  filename = paste(diff_result_dir,sample1,"_vs_",sample2,".volcano.pdf", sep='')
  cat(paste(sample1,"_vs_",sample2, sep=''),'\n')
  pdf(file=filename)
  keyvals.colour <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) >= 1,ifelse(res$log2FoldChange >= 1,'red','royalblue'),'grey')
  keyvals.colour[is.na(keyvals.colour)] <- 'grey'
  names(keyvals.colour)[keyvals.colour == 'red'] <- 'Up'
  names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'Down'
  names(keyvals.colour)[keyvals.colour == 'grey'] <- 'NS'
  
  plt <- EnhancedVolcano(res,
                         lab = rownames(res),
                         x = 'log2FoldChange',
                         y = 'pvalue',
                         colCustom = keyvals.colour,
                         xlim = c(-10, 15),
                         #ylim = c(0,30),
                         #selectLab = c('TMEM176B','ADH1A'),
                         labSize = 0, pointSize = 1,
                         # labFace = 'plain',
                         colAlpha = 0.3,
                         legendPosition = 'top', legendLabSize = 10, legendIconSize = 3.0,
                         #border = 'full',
                         borderWidth = 1.0, borderColour = 'black',
                         title = paste(sample1,"_vs_",sample2, sep=''),
                         subtitle = paste0('differential expressed genes: p<0.05, |log2FC|>1'),
                         pCutoff = 0.05,
                         FCcutoff = 1,
                         #shape = c(1, 4, 23, 25),
                         #cutoffLineType = 'solid', cutoffLineCol = 'grey', cutoffLineWidth = 2.5,
                         legendLabels=c("NS","Log2 FC","P","P & Log2 FC"),
                         #col=c('grey', 'black', 'royalblue', 'red3'),
                         )
  print(plt)
  dev.off()
}
