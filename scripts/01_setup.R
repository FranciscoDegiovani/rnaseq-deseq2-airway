install.packages("BiocManager")

library(BiocManager)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "airway", "apeglm", "pheatmap"))



# ============================================
# 1. Setup
# ============================================
library(DESeq2)
library(airway)
library(apeglm)
library(pheatmap)

# ============================================
# 2. Load data
# ============================================
data(airway)
airway
colData(airway)

# ============================================
# 3. Build DESeq2 dataset
# ============================================
dds <- DESeqDataSet(airway, design = ~ cell + dex)

# set "untrt" as the reference level (so results show trt vs untrt)
dds$dex <- relevel(dds$dex, ref = "untrt")

dds

# ============================================
# 4. Pre-filtering (remove low-count genes)
# ============================================
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds

# ============================================
# 5. Run DESeq2
# ============================================
dds <- DESeq(dds)
res <- results(dds)
summary(res)

# ============================================
# 6. Shrink log fold changes (for better visualization)
# ============================================
resultsNames(dds)  # check the exact coefficient name

res_shrunk <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm")

# MA-plot
png("results/plots/ma_plot.png", width = 800, height = 600)
plotMA(res_shrunk, ylim = c(-5, 5), main = "MA-plot: Treated vs Untreated")
dev.off()

# ============================================
# 7. Volcano plot
# ============================================
res_df <- as.data.frame(res_shrunk)
res_df$significant <- res_df$padj < 0.05 & !is.na(res_df$padj)

png("results/plots/volcano_plot.png", width = 800, height = 600)
plot(res_df$log2FoldChange, -log10(res_df$padj),
     col = ifelse(res_df$significant, "red", "grey"),
     pch = 20, cex = 0.6,
     xlab = "log2 Fold Change", ylab = "-log10 adjusted p-value",
     main = "Volcano Plot: Treated vs Untreated")
abline(v = c(-1, 1), lty = 2, col = "blue")
abline(h = -log10(0.05), lty = 2, col = "blue")
dev.off()

# ============================================
# 8. PCA plot
# ============================================
vsd <- vst(dds, blind = FALSE)

png("results/plots/pca_plot.png", width = 800, height = 600)
plotPCA(vsd, intgroup = c("dex", "cell"))
dev.off()

# ============================================
# 9. Heatmap of top significant genes
# ============================================
res_ordered <- res_df[order(res_df$padj), ]
top_genes <- rownames(res_ordered)[1:20]

mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)  # center each gene around its own mean

annotation <- as.data.frame(colData(vsd)[, c("dex", "cell")])

png("results/plots/heatmap_top20.png", width = 800, height = 900)
pheatmap(mat, annotation_col = annotation,
         main = "Top 20 Differentially Expressed Genes")
dev.off()

# ============================================
# 10. Export results table
# ============================================
res_export <- as.data.frame(res_shrunk)
res_export$gene_id <- rownames(res_export)
res_export <- res_export[order(res_export$padj), ]

write.csv(res_export, "results/deseq2_results.csv", row.names = FALSE)

sig_genes <- res_export[!is.na(res_export$padj) & res_export$padj < 0.05, ]
nrow(sig_genes)