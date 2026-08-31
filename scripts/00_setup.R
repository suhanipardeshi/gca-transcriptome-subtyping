# ==============================================================================
# Script: 00_setup.R
# Project: GCA Dissertation ML Pipeline
# Purpose: Initialize environment, set reproducibility seeds, and log versions
# Date: June 2, 2026
# Author: Suhani Balsing Pardeshi
# ==============================================================================

# 1. Reproducibility Anchor
# Enforces identical behavior for any random splitting, cross-validation, or ML modeling
set.seed(42)

# 2. Load Core Workhorse Packages
library(tidymodels)   # Master framework for modern machine learning workflows
library(ggplot2)      # Main library for generating presentation-grade figures
library(janitor)      # Used for programmatic column cleaning and quick tabulations

# 3. Document Environment for Your Dissertation Methods Section
# This automatically dumps your exact R and package versions into a text file
capture.output(sessionInfo(), file = "results/session_info_day1.txt")

cat("\n==================================================\n")
cat("SUCCESS: Environment initialized and seed set to 42.\n")
cat("Session profile exported to results/session_info_day1.txt\n")
cat("==================================================\n")

