# RNA-seq Differential Expression Analysis — DESeq2

Differential gene expression analysis using DESeq2, comparing dexamethasone-treated vs. untreated airway smooth muscle cells. Built as part of my bioinformatics/genomics portfolio, applying the same workflow used in transcriptomic research (differential expression, reproducible analysis, statistical modeling of count data).

## Data Source

- **Dataset:** airway (Bioconductor package) — RNA-seq data from Himes et al. (2014)
- **Citation:** Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, et al. (2014) RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells. PLoS ONE 9(6): e99625.
- **Design:** 4 cell lines, each with a treated and untreated sample (paired design), 8 samples total

## Tools

- R / Bioconductor
- DESeq2 (differential expression)
- apeglm (log fold-change shrinkage) — Zhu, Ibrahim & Love (2018), Bioinformatics
- pheatmap (visualization)

## Workflow

1. Load count data and sample metadata
2. Build DESeq2 dataset with design ~ cell + dex (controlling for cell line, testing treatment effect)
3. Filter low-count genes
4. Run differential expression testing
5. Shrink fold-change estimates (apeglm)
6. Visualize: MA-plot, volcano plot, PCA, heatmap of top genes
7. Export full results table (CSV)

## Results

- 4,000 genes significantly differentially expressed (adjusted p-value < 0.05) out of ~22,000 tested
- 2,610 genes up-regulated, 2,224 down-regulated with treatment
- PCA shows separation by treatment (PC2, 23% variance) alongside cell-line variability (PC1, 48% variance) — consistent with the paired experimental design

## Structure

scripts/
  analysis.R
results/
  deseq2_results.csv
  plots/
    ma_plot.png
    volcano_plot.png
    pca_plot.png
    heatmap_top20.png

## Skills Demonstrated

- RNA-seq differential expression analysis (DESeq2)
- Experimental design-aware statistical modeling (controlling for confounding variables)
- Data visualization for genomics (MA-plot, volcano plot, PCA, heatmap)
- Git version control integrated with RStudio
- Reproducible, documented analysis workflow

## Status

Complete
