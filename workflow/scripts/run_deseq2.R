library(DESeq2)
library(airway)
library(apeglm)

data(airway)
dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds$dex <- relevel(dds$dex, ref = "untrt")

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)
res_shrunk <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm")

vsd <- vst(dds, blind = FALSE)

# Save intermediate objects for the plotting step
saveRDS(dds, snakemake@output[["dds"]])
saveRDS(res_shrunk, snakemake@output[["res_shrunk"]])
saveRDS(vsd, snakemake@output[["vsd"]])

# Export results table
res_export <- as.data.frame(res_shrunk)
res_export$gene_id <- rownames(res_export)
res_export <- res_export[order(res_export$padj), ]
write.csv(res_export, snakemake@output[["csv"]], row.names = FALSE)
