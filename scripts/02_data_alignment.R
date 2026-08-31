# ==============================================================================
# Script: 02_data_alignment.R
# Purpose: Fix sample IDs, merge metadata, align samples, and save processed data
# ==============================================================================

# 1. Load the raw data (assuming 01_data_import.R code or loading directly)
source("scripts/00_setup.R")
counts_raw <- read.delim("data/raw/matrix_salmon_tximport_counts_cc_c1_c2_GENCODE.txt", row.names = 1, sep = "\t")
clinical <- read.csv("data/raw/clinical_variables_cc_c1_c2.csv")
histological <- read.csv("data/raw/histological_variables_cc_c1_c2_final.csv")

# 2. Convert the new binarized targets into factors
histological$Intima_pattern_cat  <- as.factor(histological$Intima_pattern_cat)
histological$Occlusion_grade_cat <- as.factor(histological$Occlusion_grade_cat)
histological$Media_destruction   <- as.factor(histological$Media_destruction)

# 3. Fix the "X" Prefix in RNA-seq Columns
# sub("^X", "") looks for an X at the absolute start of the string and replaces it with nothing
colnames(counts_raw) <- sub("^X", "", colnames(counts_raw))

# 4. Merge Metadata
# Combine clinical and histological data into one master dataframe using the "ID" column
meta_full <- merge(clinical, histological, by = "ID")

# 5. Identify the 79 Shared Samples
# intersect() finds only the IDs that exist in BOTH datasets
shared_samples <- intersect(colnames(counts_raw), meta_full$ID)

# 6. Subset and Force Strict Alignment
# This ensures both datasets have the exact same 79 samples in the exact same order
counts_aligned <- counts_raw[, shared_samples]

rownames(meta_full) <- meta_full$ID      # Set row names to the sample IDs
meta_aligned <- meta_full[shared_samples, ] # Order rows to match the RNA-seq columns

# 7. The Ultimate Sanity Check
# This will print TRUE if every single sample matches perfectly
is_aligned <- all(colnames(counts_aligned) == rownames(meta_aligned))
cat("\nAre the datasets perfectly aligned? :", is_aligned, "\n")

# 8. Save the Cleaned, Aligned Data for Downstream ML
# We save these to the 'processed' folder as R-specific data files (.rds)
# .rds files load incredibly fast and preserve row names and data types.
saveRDS(counts_aligned, "data/processed/counts_aligned.rds")
saveRDS(meta_aligned, "data/processed/meta_aligned.rds")

cat("Success! Aligned datasets saved to data/processed/\n")