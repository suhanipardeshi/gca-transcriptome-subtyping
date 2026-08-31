# ==============================================================================
# Script: 03_exploratory_depth.R
# Purpose: Inspect target phenotypes and check for sequencing depth imbalances
# ==============================================================================

# 1. Load setup and aligned data
source("scripts/00_setup.R")
counts <- readRDS("data/processed/counts_aligned.rds")
meta <- readRDS("data/processed/meta_aligned.rds")

# 2. Print all column names of the metadata
# We need to find the exact column names that match our 3 target phenotypes:
# (Intima_pattern_cat, Occlusion_grade_cat, Media_destruction)
cat("\n--- Metadata Columns Available ---\n")
print(colnames(meta))

# 3. Calculate Sequencing Depth per Sample
# colSums calculates the total number of read counts for each patient column
seq_depth <- colSums(counts)

# Add this sequencing depth vector to our metadata dataframe for easy plotting
meta$sequencing_depth <- seq_depth

# Summary statistics of our sequencing depths (in millions of reads)
cat("\n--- Sequencing Depth Summary (Millions of Reads) ---\n")
print(summary(meta$sequencing_depth / 1e6))

# 4. Save the updated metadata with depth column for later preprocessing
saveRDS(meta, "data/processed/meta_aligned.rds")