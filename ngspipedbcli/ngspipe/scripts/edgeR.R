library("edgeR")

# http://www.nathalievialaneix.eu/doc/html/solution_edgeR-tomato-withcode.html

data_raw <- read.table("/Users/zhangxuan/Work/Projects/bioinformatics/小桐子开花长非编码RNA/Analysis/RunAll/result/assembly_final/gene.tsv", 
                       header = TRUE, row.names="Geneid")

dim(data_raw)

head(data_raw)

cpm_log <- cpm(data_raw, log = TRUE)
median_log2_cpm <- apply(cpm_log, 1, median)
hist(median_log2_cpm)
expr_cutoff <- -1
abline(v = expr_cutoff, col = "red", lwd = 3)

sum(median_log2_cpm > expr_cutoff)

data_clean <- data_raw[median_log2_cpm > expr_cutoff, ]

cpm_log <- cpm(data_clean, log = TRUE)

heatmap(cor(cpm_log))

pca <- prcomp(t(cpm_log), scale. = TRUE)
plot(pca$x[, 1], pca$x[, 2], pch = ".", xlab = "PC1", ylab = "PC2")
text(pca$x[, 1], pca$x[, 2], labels = colnames(cpm_log))

# diff 
group <- substr(colnames(data_clean), 1, 3)
group

y <- DGEList(counts = data_clean, group = group)
y

y <- calcNormFactors(y)
y$samples

y <- estimateDisp(y)
sqrt(y$common.dispersion)
plotBCV(y)

et <- exactTest(y)
results_edgeR <- topTags(et, n = nrow(data_clean), sort.by = "none")
head(results_edgeR$table)

sum(results_edgeR$table$FDR < .1)

plotSmear(et, de.tags = rownames(results_edgeR)[results_edgeR$table$FDR < .1])
abline(h = c(-2, 2), col = "blue")



