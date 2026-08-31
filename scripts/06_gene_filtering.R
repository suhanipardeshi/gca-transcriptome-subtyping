# ==============================================================================
# Script: 06_gene_filtering.R
# Purpose: Filter low-expression and low-variance genes
# ==============================================================================

source("scripts/00_setup.R")
counts <- readRDS("data/processed/counts_aligned.rds")

cat("Original dimensions:", dim(counts)[1], "genes and", dim(counts)[2], "samples\n")

# 1. Filter Low Expression Genes
# Rule: A gene must have at least 10 counts in at least 20% of our samples (~16 patients)
min_reads <- 10
min_samples <- round(ncol(counts) * 0.20)

keep_expressed <- rowSums(counts >= min_reads) >= min_samples
counts_expressed <- counts[keep_expressed, ]

cat("Genes remaining after low-expression filter:", nrow(counts_expressed), "\n")

# 2. Filter Low Variance Genes
# Keep the top 10,000 most variable genes across all patients
gene_variances <- apply(counts_expressed, 1, var)
top_var_genes <- names(sort(gene_variances, decreasing = TRUE)[1:10000])
counts_filtered <- counts_expressed[top_var_genes, ]

cat("Genes remaining after variance filter:", nrow(counts_filtered), "\n")

# 3. Save the filtered counts
# Note: We are keeping genes in rows for now, because DESeq2 requires this format for VST normalization next!
saveRDS(counts_filtered, "data/processed/counts_filtered_top10k.rds")

cat("Success! Filtered matrix saved to data/processed/\n")