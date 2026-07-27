rule all:
    input:
        "results/deseq2_results.csv",
        "results/plots/ma_plot.png",
        "results/plots/volcano_plot.png",
        "results/plots/pca_plot.png",
        "results/plots/heatmap_top20.png"

rule run_deseq2:
    output:
        dds="results/intermediate/dds.rds",
        res_shrunk="results/intermediate/res_shrunk.rds",
        vsd="results/intermediate/vsd.rds",
        csv="results/deseq2_results.csv"
    script:
        "workflow/scripts/run_deseq2.R"

rule make_plots:
    input:
        dds="results/intermediate/dds.rds",
        res_shrunk="results/intermediate/res_shrunk.rds",
        vsd="results/intermediate/vsd.rds"
    output:
        ma="results/plots/ma_plot.png",
        volcano="results/plots/volcano_plot.png",
        pca="results/plots/pca_plot.png",
        heatmap="results/plots/heatmap_top20.png"
    script:
        "workflow/scripts/make_plots.R"