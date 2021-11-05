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
  'k_anno'   , 't', 1, "character", 'k number from kaas or eggnog',
  'output'  , 'o', 1, "character", 'enriched table file'
), byrow=TRUE, ncol=5)
args = getopt(command)

if ( !is.null(args$help) || is.null(args$k_anno) || is.null(args$output) ) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}

chooseCRANmirror(ind=16)
# install normal packages
my_packages_normal <- c("devtools", "ggplot2", "BiocManager")                                        # Specify your packages
not_installed_normal <- my_packages_normal[!(my_packages_normal %in% installed.packages()[ , "Package"])]    # Extract not installed packages
if(length(not_installed_normal)) install.packages(not_installed_normal)                               # Install not installed packages

# install bioconductor packages
my_packages_bioconductor <- c("clusterProfiler")                                        # Specify your packages
not_installed_bioconductor <- my_packages_bioconductor[!(my_packages_bioconductor %in% installed.packages()[ , "Package"])]    # Extract not installed packages
if(length(not_installed_bioconductor)) install.packages(not_installed_bioconductor)                               # Install not installed packages

#devtools::install_github('GuangchuangYu/clusterProfiler')

library(clusterProfiler)


# clusterProfiler deal with go
# go2term('GO:0043232')
# go2ont('GO:0043232')

# clusterProfiler deal with kegg
# bitr_kegg(c("K05692","K19036"), "kegg", "Path", "ko") -> x
# ko2name(x$Path) -> y
# merge(x, y, by.x='Path', by.y='ko')

k <- read.table(args$k_anno, sep="\t", quote="")

bitr_kegg(k$V2, "kegg", "Path", "ko") -> x

unique(ko2name(x$Path)) -> y

kegg <- merge(x, y, by.x='Path', by.y='ko')

anno_finnal <- unique(merge(k, kegg, by.k='V2', by.kegg='kegg'))

write.table(unique(anno_finnal[c('V1', 'Path', 'name')]), file=args$output, col.names = FALSE, sep='\t', quote = FALSE, row.names = FALSE)



