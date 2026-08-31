# ==============================================================================
# Script: 23_figures.R
# Purpose: Generate dissertation figures from the benchmark and diagnostic runs.
# ==============================================================================
setwd("C:/Users/suhan/Desktop/Dissertation_aug")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })

FIG <- "figures/final"; dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
RES <- "results/benchmark"

ph_lab <- c(intima = "Intimal inflammation", occlusion = "Occlusion grade",
            media = "Media destruction")
al_lab <- c(lasso = "LASSO", ridge = "Ridge", enet = "Elastic net",
            svm = "SVM (RBF)", rf = "Random forest", xgb = "XGBoost")

theme_set(theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom"))

# ---------------------------------------------------------------- Figure 1
folds <- do.call(rbind, lapply(list.files(RES, "^folds_.*csv$", full.names = TRUE),
  function(f) {
    p <- sub("^folds_([a-z]+)_([a-z]+)\\.csv$", "\\1|\\2", basename(f))
    d <- read.csv(f)
    d$phenotype <- strsplit(p, "\\|")[[1]][1]
    d$algorithm <- strsplit(p, "\\|")[[1]][2]
    d
  }))
folds$phenotype <- factor(ph_lab[folds$phenotype], levels = ph_lab)
folds$algorithm <- factor(al_lab[folds$algorithm], levels = al_lab)

smry <- folds %>% group_by(phenotype, algorithm) %>%
  summarise(m = mean(roc_auc, na.rm = TRUE), s = sd(roc_auc, na.rm = TRUE),
            .groups = "drop")

p1 <- ggplot(smry, aes(algorithm, m, fill = algorithm)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40") +
  geom_col(width = .68, alpha = .9) +
  geom_errorbar(aes(ymin = pmax(0, m - s), ymax = pmin(1, m + s)), width = .18) +
  geom_jitter(data = folds, aes(algorithm, roc_auc), width = .13,
              size = .5, alpha = .35, inherit.aes = FALSE) +
  facet_wrap(~phenotype) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "ROC-AUC (nested CV)") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
ggsave(file.path(FIG, "fig1_benchmark_auc.pdf"), p1, width = 9, height = 4)
ggsave(file.path(FIG, "fig1_benchmark_auc.png"), p1, width = 9, height = 4, dpi = 160)

# ---------------------------------------------------------------- Figure 2
counts <- readRDS("processed/counts_aligned.rds")
vst    <- readRDS("processed/vst_counts_ml.rds")
meta   <- readRDS("processed/meta_aligned.rds")
meta$cohort <- sub(".*_", "", meta$ID)
meta$cohort <- factor(meta$cohort, levels = c("c1", "c2", "cc"),
                      labels = c("Cohort 1", "Cohort 2", "Concatenated"))
meta$depth <- colSums(counts) / 1e6

pca <- prcomp(vst, scale. = FALSE)
ve  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pc  <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                  cohort = meta$cohort, depth = meta$depth)

p2a <- ggplot(pc, aes(PC1, PC2, colour = cohort)) +
  geom_point(size = 2, alpha = .85) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
       colour = NULL, title = "A. Batch structure")

p2b <- ggplot(pc, aes(depth, PC2, colour = cohort)) +
  geom_point(size = 2, alpha = .85) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey25",
              linewidth = .5, inherit.aes = FALSE, aes(depth, PC2)) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Sequencing depth (millions of reads)",
       y = sprintf("PC2 (%.1f%%)", ve[2]), colour = NULL,
       title = "B. PC2 tracks sequencing depth")

dd <- data.frame(depth = meta$depth,
                 media = factor(as.character(meta$Media_destruction),
                                levels = c("0", "1"),
                                labels = c("Absent", "Present")))
dd <- dd[!is.na(dd$media), ]
p2c <- ggplot(dd, aes(media, depth, fill = media)) +
  geom_boxplot(width = .55, alpha = .8, outlier.shape = NA) +
  geom_jitter(width = .14, size = 1, alpha = .6) +
  scale_fill_brewer(palette = "Set1", guide = "none") +
  labs(x = "Media destruction", y = "Sequencing depth (millions)",
       title = "C. Depth differs by phenotype")

pdf(file.path(FIG, "fig2_confounding.pdf"), width = 11, height = 3.6)
gridExtra::grid.arrange(p2a, p2b, p2c, ncol = 3)
dev.off()
png(file.path(FIG, "fig2_confounding.png"), width = 11, height = 3.6,
    units = "in", res = 160)
gridExtra::grid.arrange(p2a, p2b, p2c, ncol = 3)
dev.off()

# ---------------------------------------------------------------- Figure 3
roc_pts <- function(prob, y) {
  o <- order(prob, decreasing = TRUE); y <- y[o]
  tpr <- cumsum(y) / sum(y); fpr <- cumsum(1 - y) / sum(1 - y)
  data.frame(fpr = c(0, fpr), tpr = c(0, tpr))
}
best <- c(intima = "ridge", occlusion = "rf", media = "ridge")
rocs <- do.call(rbind, lapply(names(best), function(p) {
  d <- read.csv(file.path(RES, sprintf("preds_%s_%s.csv", p, best[p])))
  r <- roc_pts(d$prob, d$truth)
  r$phenotype <- ph_lab[p]
  r$lab <- sprintf("%s (%s)", ph_lab[p], al_lab[best[p]])
  r
}))
p3 <- ggplot(rocs, aes(fpr, tpr, colour = lab)) +
  geom_abline(linetype = 2, colour = "grey55") +
  geom_line(linewidth = .8) +
  scale_colour_brewer(palette = "Dark2") +
  coord_equal() +
  labs(x = "False positive rate", y = "True positive rate", colour = NULL)
ggsave(file.path(FIG, "fig3_roc_best.pdf"), p3, width = 5.5, height = 5.2)
ggsave(file.path(FIG, "fig3_roc_best.png"), p3, width = 5.5, height = 5.2, dpi = 160)

cat("Figures 1-3 written to", FIG, "\n")

# ---------------------------------------------------------------- Figure 4
pf <- file.path(RES, "permutation_distribution.csv")
if (file.exists(pf)) {
  perm <- read.csv(pf)
  obs  <- read.csv(file.path(RES, "baselines_and_nulls.csv"))
  perm$phenotype <- factor(ph_lab[perm$phenotype], levels = ph_lab)
  obs$phenotype  <- factor(ph_lab[obs$phenotype],  levels = ph_lab)
  p4 <- ggplot(perm, aes(perm_auc)) +
    geom_histogram(bins = 30, fill = "grey70", colour = "white") +
    geom_vline(data = obs, aes(xintercept = lasso_auc), colour = "firebrick",
               linewidth = .9) +
    facet_wrap(~phenotype, scales = "free_y") +
    labs(x = "ROC-AUC under permuted labels", y = "Frequency")
  ggsave(file.path(FIG, "fig4_permutation.pdf"), p4, width = 9, height = 3.2)
  ggsave(file.path(FIG, "fig4_permutation.png"), p4, width = 9, height = 3.2, dpi = 160)
  cat("Figure 4 written.\n")
} else cat("Permutation results not yet available - rerun for Figure 4.\n")
