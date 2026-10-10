# Install the packages used by "RNAseqanalysisin BC.R".
# Run once with:  Rscript install.R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db"), update = FALSE, ask = FALSE)
install.packages(c("glmnet", "pheatmap", "survival"))
