library("ggplot2")
library("pheatmap")
library("reshape2")
library("stringr")

# read expression matrix
data <-read.table("../result/timeSerise/IB1_IB2_IB3_FFB_FF.matrix", header = T,row.names=1)
head(data,10L)

# read tissueEnriched gene
ts <- read.table("../result/tissueEnrich/ts.tsv", header = T)

#ts_matrix = subset(data, rownames(data) %in% unique(ts$Gene) & str_detect(rownames(data), "XLOC") ) 
ts_matrix = subset(data, rownames(data) %in% unique(ts$Gene)) 

annotation_row = colnames(data)

# no cluster number set up
p<-pheatmap(ts_matrix,show_rownames=FALSE,scale = "column",clustering_distance_rows = "correlation", 
            color = colorRampPalette(c("navy", "white", "firebrick3"))(50),cluster_cols = FALSE,cluster_rows=T,legend = FALSE)

# set up cluster numbers
p<-pheatmap(ts_matrix,cutree_rows =10, show_rownames=FALSE,scale = "column", 
            color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
            cluster_cols = FALSE,cluster_rows=T,legend = FALSE)

# save cluster genes
list = pheatmap(ts_matrix, cutree_rows =10, scale = "column",show_rownames=FALSE, 
                color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                cluster_cols = FALSE,cluster_rows=T,legend = FALSE)

row_cluster = cutree(list$tree_row,k=10)
head(row_cluster)
newOrder=data[list$tree_row$order,]
newOrder[,ncol(newOrder)+1]=row_cluster[match(rownames(newOrder),names(row_cluster))]
colnames(newOrder)[ncol(newOrder)]="Cluster"
head(newOrder)
write.csv(newOrder,"XXX.csv")

