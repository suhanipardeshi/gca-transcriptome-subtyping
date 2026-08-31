# ==============================================================================
# Script: 29_gsea_kegg.R
# Purpose: Complete the functional analysis with (a) KEGG over-representation,
#          and (b) GSEA on the full ranked gene list.
#
# Why GSEA rather than only over-representation (ORA):
#   ORA takes a SELECTED gene list and asks whether any pathway is over-
#   represented in it. That is the wrong tool here for two reasons established
#   earlier in this project: the LASSO signatures are unstable (median gene
#   retained in 1 of 15 folds), and random gene sets of the same size predict
#   as well as the fitted ones. ORA applied to such a list is testing an
#   arbitrary draw.
#
#   GSEA instead ranks ALL 10,000 genes by their association with the phenotype
#   and asks whether members of a pathway cluster towards either extreme. It
#   requires no gene selection at all, so it is the natural test of the central
#   finding that the signal is diffuse: if many genes each contribute a little,
#   a coordinated pathway shift should be detectable even though no individual
#   gene survives selection.
#
# Ranking statistic: moderated t-statistic (limma) for severe vs mild, computed
# on the VST matrix. This is the univariate association of each gene with the
# phenotype, which is also the quantity the maximal-shrinkage ridge model
# weights genes by - so the ranking is aligned with the best-performing model.
# ==============================================================================

setwd("C:/Users/suhan/Desktop/Dissertation_aug")
suppressPackageStartupMessages({
  library(limma); library(clusterProfiler); library(org.Hs.eg.db)
  library(AnnotationDbi)
})
RES <- "results/benchmark"; dir.create(RES, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

vst  <- readRDS("processed/vst_counts_ml.rds")
meta <- readRDS("processed/meta_aligned.rds")
num01 <- function(x) as.numeric(as.character(x))
universe_ens <- colnames(vst)

PH <- list(
  intima    = list(lab = "Intima pattern",    y = num01(meta$Intima_pattern_cat),
                   cc = !is.na(meta$Intima_pattern)),
  occlusion = list(lab = "Occlusion grade",   y = num01(meta$Occlusion_grade_cat),
                   cc = !is.na(meta$Occlusion_grade)),
  media     = list(lab = "Media destruction", y = num01(meta$Media_destruction),
                   cc = !is.na(num01(meta$Media_destruction))))

# ENSEMBL -> ENTREZ map for the whole modelled universe (KEGG needs Entrez)
map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = universe_ens,
         keytype = "ENSEMBL", columns = "ENTREZID"))
map <- map[!is.na(map$ENTREZID) & !duplicated(map$ENSEMBL), ]
cat("universe: ", length(universe_ens), " ENSEMBL genes; ",
    nrow(map), " map to ENTREZ\n\n", sep = "")

summary_rows <- list()

for (nm in names(PH)) {
  ph <- PH[[nm]]
  keep <- ph$cc & !is.na(ph$y)
  X <- t(as.matrix(vst[keep, , drop = FALSE]))   # limma wants genes x samples
  y <- factor(ifelse(ph$y[keep] == 1, "severe", "mild"), levels = c("mild", "severe"))

  cat("==========================================================\n")
  cat(ph$lab, " (n = ", ncol(X), "; severe = ", sum(y == "severe"), ")\n", sep = "")
  cat("==========================================================\n")

  # ---- rank all genes by moderated t-statistic ---------------------------
  design <- model.matrix(~ y)
  fit <- eBayes(lmFit(X, design))
  tt  <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  cat("genes with FDR < 0.05 (univariate):", sum(tt$adj.P.Val < 0.05), "\n")

  rank_ens <- sort(setNames(tt$t, rownames(tt)), decreasing = TRUE)
  write.csv(data.frame(gene = names(rank_ens), t = as.numeric(rank_ens)),
            file.path(RES, sprintf("ranked_genes_%s.csv", nm)), row.names = FALSE)

  # ---- GSEA: GO biological process ---------------------------------------
  g <- tryCatch(gseGO(geneList = rank_ens, OrgDb = org.Hs.eg.db,
                      keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.05,
                      pAdjustMethod = "BH", minGSSize = 15, maxGSSize = 500,
                      seed = TRUE, verbose = FALSE),
                error = function(e) { cat("  gseGO error:", conditionMessage(e), "\n"); NULL })
  ngo <- if (is.null(g)) 0 else nrow(as.data.frame(g))
  cat("GSEA GO BP  : ", ngo, " terms at FDR < 0.05\n", sep = "")
  if (ngo > 0) {
    d <- as.data.frame(g)
    d <- d[order(d$p.adjust), ]
    write.csv(d, file.path(RES, sprintf("gsea_go_%s.csv", nm)), row.names = FALSE)
    print(head(d[, c("Description", "setSize", "NES", "p.adjust")], 10), row.names = FALSE)
  }

  # ---- GSEA: KEGG ---------------------------------------------------------
  rk <- rank_ens[names(rank_ens) %in% map$ENSEMBL]
  names(rk) <- map$ENTREZID[match(names(rk), map$ENSEMBL)]
  rk <- sort(rk[!duplicated(names(rk))], decreasing = TRUE)
  k <- tryCatch(gseKEGG(geneList = rk, organism = "hsa", pvalueCutoff = 0.05,
                        pAdjustMethod = "BH", minGSSize = 15, maxGSSize = 500,
                        seed = TRUE, verbose = FALSE),
                error = function(e) { cat("  gseKEGG error:", conditionMessage(e), "\n"); NULL })
  nkg <- if (is.null(k)) 0 else nrow(as.data.frame(k))
  cat("GSEA KEGG   : ", nkg, " pathways at FDR < 0.05\n", sep = "")
  if (nkg > 0) {
    d <- as.data.frame(k); d <- d[order(d$p.adjust), ]
    write.csv(d, file.path(RES, sprintf("gsea_kegg_%s.csv", nm)), row.names = FALSE)
    print(head(d[, c("Description", "setSize", "NES", "p.adjust")], 10), row.names = FALSE)
  }

  # ---- KEGG over-representation on the selected genes (for completeness) --
  sf <- file.path(RES, sprintf("stability_%s.csv", nm))
  nora <- NA
  if (file.exists(sf)) {
    sel <- read.csv(sf)$gene
    sel_e <- map$ENTREZID[match(sel, map$ENSEMBL)]
    sel_e <- sel_e[!is.na(sel_e)]
    ko <- tryCatch(enrichKEGG(gene = sel_e, universe = map$ENTREZID,
                              organism = "hsa", pvalueCutoff = 0.05,
                              pAdjustMethod = "BH"),
                   error = function(e) NULL)
    nora <- if (is.null(ko)) 0 else nrow(as.data.frame(ko))
    cat("KEGG ORA on ", length(sel_e), " ever-selected genes: ", nora,
        " pathways at FDR < 0.05\n", sep = "")
    if (nora > 0) {
      d <- as.data.frame(ko)
      write.csv(d, file.path(RES, sprintf("kegg_ora_%s.csv", nm)), row.names = FALSE)
      print(head(d[, c("Description", "GeneRatio", "p.adjust")], 8), row.names = FALSE)
    }
  }
  cat("\n")
  summary_rows[[nm]] <- data.frame(phenotype = nm, label = ph$lab,
    n = ncol(X), univariate_fdr05 = sum(tt$adj.P.Val < 0.05),
    gsea_go_terms = ngo, gsea_kegg_paths = nkg, kegg_ora_paths = nora)
}

s <- do.call(rbind, summary_rows)
write.csv(s, file.path(RES, "functional_analysis_summary.csv"), row.names = FALSE)
cat("================ FUNCTIONAL ANALYSIS SUMMARY ================\n")
print(s, row.names = FALSE)
