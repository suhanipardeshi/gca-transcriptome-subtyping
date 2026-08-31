# ==============================================================================
# Script: 25_enrichment_universe.R
# Purpose: (a) annotate the most stably selected genes;
#          (b) quantify how much the choice of background universe inflates
#              enrichment significance.
#
# The original pipeline called enrichGO() without a `universe` argument, so the
# background defaulted to every annotated human gene rather than the 10,000
# variance-filtered genes actually eligible for selection. This script runs the
# same gene lists both ways so the difference can be reported, not assumed.
# ==============================================================================
setwd("C:/Users/suhan/Desktop/Dissertation_aug")
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)
})
RES <- "results/benchmark"
counts <- readRDS("processed/vst_counts_ml.rds")
universe <- colnames(counts)
cat("background universe (genes entering models):", length(universe), "\n\n")

PH <- c(intima = "Intimal inflammation severity",
        occlusion = "Arterial occlusion grade",
        media = "Media destruction")

rows <- list()
for (nm in names(PH)) {
  s <- read.csv(file.path(RES, sprintf("stability_%s.csv", nm)))
  s$symbol <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = s$gene,
                column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  write.csv(s, file.path(RES, sprintf("stability_%s_annotated.csv", nm)),
            row.names = FALSE)

  cat("=====", PH[nm], "=====\n")
  cat("Most stably selected genes (of 15 outer folds):\n")
  print(head(s[, c("symbol", "gene", "n_folds", "frequency")], 10),
        row.names = FALSE)

  gl <- s$gene                                  # union of ever-selected genes
  run_go <- function(univ, tag) {
    go <- tryCatch(enrichGO(gene = gl, universe = univ, OrgDb = org.Hs.eg.db,
             keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH",
             pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
             error = function(e) NULL)
    d <- if (is.null(go)) data.frame() else as.data.frame(go)
    cat(sprintf("  %-28s significant GO BP terms: %d\n", tag, nrow(d)))
    d
  }
  cat("\nGO over-representation on", length(gl), "ever-selected genes:\n")
  d_correct <- run_go(universe, "correct universe (10,000)")
  d_wrong   <- run_go(NULL,     "default whole-genome universe")

  if (nrow(d_correct)) {
    write.csv(d_correct, file.path(RES, sprintf("go_correct_%s.csv", nm)),
              row.names = FALSE)
    print(head(d_correct[, c("Description", "GeneRatio", "pvalue", "p.adjust")], 5),
          row.names = FALSE)
  }
  if (nrow(d_wrong)) {
    write.csv(d_wrong, file.path(RES, sprintf("go_wronguniverse_%s.csv", nm)),
              row.names = FALSE)
    cat("  (top term under the inflated background: ",
        d_wrong$Description[1], ", p.adj = ",
        signif(d_wrong$p.adjust[1], 3), ")\n", sep = "")
  }
  rows[[nm]] <- data.frame(phenotype = nm, n_genes = length(gl),
                           n_terms_correct = nrow(d_correct),
                           n_terms_wrong = nrow(d_wrong))
  cat("\n")
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path(RES, "enrichment_universe_comparison.csv"), row.names = FALSE)
cat("===== UNIVERSE COMPARISON =====\n"); print(tab, row.names = FALSE)
