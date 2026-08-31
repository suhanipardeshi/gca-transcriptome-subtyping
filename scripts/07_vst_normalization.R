# ==============================================================================
# Script: 07_vst_normalization.R
# Purpose: Apply DESeq2 VST to remove depth bias and format for Machine Learning
# ==============================================================================

source("scripts/00_setup.R")
library(DESeq2)

counts <- readRDS("data/processed/counts_filtered_top10k.rds")
meta <- readRDS("data/processed/meta_aligned.rds")

# FIX: Round Salmon's estimated decimal counts to whole integers for DESeq2
counts <- round(counts)

# 1. Final Safety Check
stopifnot(all(colnames(counts) == meta$ID))

# 2. Build the DESeq2 Dataset Object
cat("\nBuilding DESeqDataSet...\n")
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~ 1)

# 3. Apply the Variance Stabilizing Transformation (VST)
cat("Running Variance Stabilizing Transformation (VST). This may take a minute...\n")
vsd <- vst(dds, blind = TRUE)

# 4. Extract the normalized matrix
vst_matrix <- assay(vsd)

# 5. Transpose for Machine Learning (Rows = Patients, Columns = Genes)
counts_ml_ready <- as.data.frame(t(vst_matrix))

# 6. Save the ultimate ML-ready dataset
saveRDS(counts_ml_ready, "data/processed/vst_counts_ml.rds")

cat("\nSUCCESS: VST normalization complete.\n")
cat("Final ML Matrix Dimensions:", dim(counts_ml_ready)[1], "patients x", dim(counts_ml_ready)[2], "genes\n")
