# GCA Transcriptome Subtyping

R code for an MSc dissertation on predicting histological phenotypes of giant cell arteritis (GCA) from temporal artery biopsy RNA sequencing data. The dissertation is titled "Robust ML Driven Subtyping of Giant Cell Arteritis from Temporal Artery Transcriptomes."

## What this repository contains

Only the scripts. No raw or processed patient data, and no per sample results, are included here. The source data is de identified clinical and histological information from the UK GCA Consortium, and analysis was confined to a secure institutional computing environment as stated in the dissertation methods. This repository exists to make the analysis code available, not the data.

## Pipeline overview

The scripts run in the following order.

1. `00_setup.R`. Loads required packages and sets a fixed seed for reproducibility.
2. `01_data_import.R`. Imports raw RNA sequencing counts and clinical metadata.
3. `02_data_alignment.R`. Aligns sample identifiers between the count matrix and metadata, and saves the aligned data.
4. `03_exploratory_depth.R`. Computes sequencing depth per sample and adds it to the metadata. This step is required by later confounding checks.
5. `06_gene_filtering.R`. Filters low expression and low variance genes down to the 10,000 most variable genes.
6. `07_vst_normalization.R`. Applies a DESeq2 variance stabilising transformation to the filtered counts.
7. `21_benchmark_fast.R`. Benchmarks six algorithms (LASSO, ridge, elastic net, support vector machine, random forest, XGBoost) for three phenotypes (intima pattern, occlusion grade, media destruction) under repeated stratified nested cross validation.
8. `22_baselines_nulls.R`. Computes non transcriptomic baselines, depth adjusted models, random gene set nulls, and label permutation nulls for LASSO.
9. `23_figures.R`. Generates the benchmark performance figure, the confounding figure, and the ROC curve figure.
10. `24_stability_enrichment.R`. Refits LASSO across all outer folds to measure how often each gene is selected, and runs over representation analysis on the resulting gene lists.
11. `25_enrichment_universe.R` and `26_universe_test_original.R`. Check how the choice of background gene universe affects over representation results, including a reproducibility check against an earlier stored signature.
12. `27_ridge_limit_permutation.R`. Evaluates ridge regression at fixed penalty values across ten orders of magnitude, and runs a label permutation test on ridge specifically.
13. `28_figures_nulls.R`. Generates the permutation null distribution figure.
14. `29_gsea_kegg.R`. Runs gene set enrichment analysis against Gene Ontology and KEGG using the full ranked gene list, with no gene selection step.
15. `30_figures_gsea.R`. Generates the gene set enrichment analysis figures.
16. `31_feature_counts.R`. Compares how many genes each of the six algorithms actually uses, as a supplementary sparsity comparison across algorithms. This uses a single fixed hyperparameter per phenotype rather than the fold by fold retuning used in script 24, so its LASSO gene counts differ from script 24 by design and are not directly comparable to it.

## Requirements

R 4.5.1 or later, with the following packages: glmnet, ranger, xgboost, kernlab, DESeq2, clusterProfiler, org.Hs.eg.db.

## Note on excluded scripts

The dissertation project originally included an earlier pipeline built on a single train test split rather than nested cross validation. That earlier pipeline is not included here because it was superseded and its results are not the ones reported in the dissertation.
