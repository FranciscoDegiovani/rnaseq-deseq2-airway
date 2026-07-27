library(DESeq2)
library(pheatmap)

dds <- readRDS(snakemake@input[["dds"]])
res_shrunk <- readRDS(snakemake@input[["res_shrunk"]])
vsd <- readRDS(snakemake@input[["vsd"]])

# MA-plot
png(snakemake@output[["ma"]], width = 800, height = 600)
plotMA(res_shrunk, ylim = c(-5, 5), main = "MA-plot: Treated vs Untreated")
dev.off()

# Volcano plot
res_df <- as.data.frame(res_shrunk)
res_df$significant <- res_df$padj < 0.05 & !is.na(res_df$padj)

png(snakemake@output[["volcano"]], width = 800, height = 600)
plot(res_df$log2FoldChange, -log10(res_df$padj),
     col = ifelse(res_df$significant, "red", "grey"),
     pch = 20, cex = 0.6,
     xlab = "log2 Fold Change", ylab = "-log10 adjusted p-value",
     main = "Volcano Plot: Treated vs Untreated")
abline(v = c(-1, 1), lty = 2, col = "blue")
abline(h = -log10(0.05), lty = 2, col = "blue")
dev.off()

# PCA
png(snakemake@output[["pca"]], width = 800, height = 600)
plotPCA(vsd, intgroup = c("dex", "cell"))
dev.off()

# Heatmap of top 20 genes
res_ordered <- res_df[order(res_df$padj), ]
top_genes <- rownames(res_ordered)[1:20]
mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)
annotation <- as.data.frame(colData(vsd)[, c("dex", "cell")])

png(snakemake@output[["heatmap"]], width = 800, height = 900)
pheatmap(mat, annotation_col = annotation, main = "Top 20 Differentially Expressed Genes")
dev.off()
