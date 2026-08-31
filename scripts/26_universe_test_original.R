# ==============================================================================
# Script: 26_universe_test_original.R
# Purpose: Direct test of the background-universe effect on the ORIGINAL
#          single-split signatures (those in 'Revised results/'), which are the
#          lists that produced the previously reported enrichment terms.
#          Each list is analysed twice: once against the whole annotated genome
#          (the clusterProfiler default, as originally used) and once against
#          the 10,000 variance-filtered genes actually eligible for selection.
# ==============================================================================
setwd("C:/Users/suhan/Desktop/Dissertation_aug")
suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db)})

universe <- colnames(readRDS("processed/vst_counts_ml.rds"))
files <- c(
  "LASSO intima"    = "Revised results/lasso_intima_34_genes.csv",
  "LASSO occlusion" = "Revised results/lasso_occlusion_genes.csv",
  "LASSO media"     = "Revised results/lasso_media_genes.csv",
  "ENet intima"     = "Revised results/enet_intima_genes.csv",
  "ENet occlusion"  = "Revised results/enet_occlusion_genes.csv",
  "ENet media"      = "Revised results/enet_media_genes.csv")

go_n <- function(gl, univ) {
  r <- tryCatch(enrichGO(gene = gl, universe = univ, OrgDb = org.Hs.eg.db,
        keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH",
        pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
        error = function(e) NULL)
  if (is.null(r)) return(data.frame())
  as.data.frame(r)
}

rows <- list()
for (nm in names(files)) {
  gl <- read.csv(files[nm], stringsAsFactors = FALSE)$Variable
  dw <- go_n(gl, NULL)        # default: whole annotated genome (as originally run)
  dc <- go_n(gl, universe)    # correct: the 10,000 genes eligible for selection
  cat(sprintf("%-16s n=%3d | whole-genome background: %3d terms | correct background: %3d terms\n",
              nm, length(gl), nrow(dw), nrow(dc)))
  if (nrow(dw)) cat("      top (whole-genome): ", dw$Description[1],
                    "  p.adj=", signif(dw$p.adjust[1], 3), "\n", sep = "")
  if (nrow(dc)) cat("      top (correct)     : ", dc$Description[1],
                    "  p.adj=", signif(dc$p.adjust[1], 3), "\n", sep = "")
  rows[[nm]] <- data.frame(signature = nm, n_genes = length(gl),
                           terms_wholegenome = nrow(dw), terms_correct = nrow(dc))
}
tab <- do.call(rbind, rows)
write.csv(tab, "results/benchmark/universe_test_original.csv", row.names = FALSE)
cat("\n===== SUMMARY =====\n"); print(tab, row.names = FALSE)
