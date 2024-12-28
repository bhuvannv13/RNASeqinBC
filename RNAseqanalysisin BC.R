
  # Load necessary libraries

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

if (!require("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")

if (!require("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")

if (!require("pheatmap", quietly = TRUE))
  install.packages("pheatmap")

if (!require("glmnet", quietly = TRUE))
  install.packages("glmnet")

library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(glmnet)


## Load and Preprocess Data

### Set File Paths

# Define the paths for the datasets
folder_path <- "C:/Users/bhuva/Downloads/brca_tcga_pan_can_atlas_2018"  # Replace with actual extraction path
data_clinical_path <- paste(folder_path, "data_clinical_patient.txt", sep = "/")
data_rna_path <- paste(folder_path, "data_mrna_seq_v2_rsem.txt", sep = "/")
data_cna_path <- paste(folder_path, "data_cna.txt", sep = "/")
head(data_cna)
head(data_rna_seq)
head(data_patient)
### Read Patient Clinical Data

clinical_data <- read.delim(data_clinical_path, header = TRUE, sep = "\t")
# Skip first rows if they are descriptions
clinical_data <- clinical_data[5:nrow(clinical_data), ]

### Read RNA-Seq Data

rna_data <- read.delim(data_rna_path, header = TRUE, sep = "\t")
# Extract gene names and ensure matrix format
assay <- round(as.matrix(rna_data[, -c(1, 2)]))
rownames(assay) <- rna_data[, 1]


### Read CNA Data

cna_data <- read.delim(data_cna_path, header = TRUE, sep = "\t")

#### Analyse the data
unique(data_cna)
unique(data_rna_seq)

colnames(data_cna)

# Verify dimensions
print(dim(clinical_data))

ERBB2_row = which(data_cna[, 1] == "ERBB2")
ERBB2_data = matrix(data_cna[ERBB2_row, 3:ncol(data_cna)], 
                    ncol = 1, 
                    dimnames = list(colnames(data_cna)[3:ncol(data_cna)], "ERBB2_Count"))

## Match Patient IDs and Create Metadata

# Extract patient barcodes and map to clinical data
pat_ids <- clinical_data[, 1]
metadata <- matrix(0, ncol = 2, nrow = ncol(assay))
colnames(metadata) <- c("HER2_Amplified", "Stage")

for (i in 1:ncol(assay)) {
  pat_barcode <- substr(colnames(assay)[i], 1, 12)
  pat_barcode <- gsub("\\.", "-", pat_barcode)
  idx <- which(pat_ids == pat_barcode)
  metadata[i, 1] <- as.numeric(cna_data[idx, "ERBB2"] > 0)  # HER2 Amplification
  metadata[i, 2] <- clinical_data[idx, "Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code"]
}


## Normalize Data and Differential Expression Analysis

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = assay,
                              colData = as.data.frame(metadata),
                              design = ~ HER2_Amplified)

# Pre-filter low-count genes
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]

# Normalize and run differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

# Extract top 10 differentially expressed genes
res <- res[order(res$padj), ]
top_genes <- head(res, 10)


## Pathway Enrichment Analysis

# Subset significant genes
significant_genes <- rownames(res[res$padj < 0.05, ])

# Separate over and under-expressed genes
over_expressed <- significant_genes[res$log2FoldChange > 0]
under_expressed <- significant_genes[res$log2FoldChange < 0]

# Perform GO enrichment
go_results <- enrichGO(
  gene = over_expressed,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH"
)

# Plot enrichment results
dotplot(go_results, showCategory = 10) + ggtitle("GO Enrichment for Over-Expressed Genes")


## PCA and Heatmap

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# PCA plot
plotPCA(vsd, intgroup = c("HER2_Amplified"))

# Heatmap for top differentially expressed genes
vsd_top <- assay(vsd)[rownames(top_genes), ]
pheatmap(vsd_top, cluster_rows = TRUE, cluster_cols = TRUE, scale = 'row')


## Survival Model with Lasso Regularization

# Prepare data for Cox regression
survival_data <- clinical_data[, c("Overall.Survival..Months.", "Overall.Survival.Status")]
expression_data <- assay(vsd)[rownames(top_genes), ]

# Fit Lasso regression
cv_fit <- cv.glmnet(x = t(expression_data), y = survival_data[, 1], family = "cox")
plot(cv_fit)
