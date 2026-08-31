# ==============================================================================
# Script: 01_data_import.R
# Purpose: Import GCA raw counts and metadata, perform initial sanity checks
# ==============================================================================

# 1. Load setup script 
source("scripts/00_setup.R")

# 2. Import the Raw Gene Expression Counts
counts_raw <- read.delim("data/raw/matrix_salmon_tximport_counts_cc_c1_c2_GENCODE.txt", 
                         header = TRUE, 
                         row.names = 1,
                         sep = "\t")

# 3. Import the Clinical and Histological Metadata
# Using base R read.csv since they are standard comma-separated files
clinical <- read.csv("data/raw/clinical_variables_cc_c1_c2.csv")
histological <- read.csv("data/raw/histological_variables_cc_c1_c2_final.csv")

# 4. Initial Sanity Checks (Dimensions)
cat("\n--- Data Dimensions ---\n")
cat("Counts (Genes x Samples):", dim(counts_raw), "\n")
cat("Clinical (Samples x Variables):", dim(clinical), "\n")
cat("Histological (Samples x Variables):", dim(histological), "\n")

# Convert the newly binarized targets into factors
histological$Intima_pattern_cat <- as.factor(histological$Intima_pattern_cat)
histological$Occlusion_grade_cat <- as.factor(histological$Occlusion_grade_cat)
histological$Media_destruction <- as.factor(histological$Media_destruction)