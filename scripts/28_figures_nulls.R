# ==============================================================================
# Script: 28_figures_nulls.R
# Purpose: Figures for the permutation nulls and the ridge penalty plateau.
# ==============================================================================
setwd("C:/Users/suhan/Desktop/Dissertation_aug")
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
FIG <- "figures/final"; RES <- "results/benchmark"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

ph_lab <- c(intima = "Intima pattern", occlusion = "Occlusion grade",
            media = "Media destruction")
theme_set(theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA)))

# ---- Figure 4: permutation null vs observed (ridge) -----------------------
perm <- read.csv(file.path(RES, "ridge_permutation_distribution.csv"))
obs  <- read.csv(file.path(RES, "ridge_permutation.csv"))
perm$phenotype <- factor(ph_lab[perm$phenotype], levels = ph_lab)
obs$phenotype  <- factor(ph_lab[obs$phenotype],  levels = ph_lab)
obs$lab <- sprintf("observed %.3f\np = %.3f", obs$observed_auc, obs$perm_p)

p4 <- ggplot(perm, aes(perm_auc)) +
  geom_histogram(bins = 28, fill = "grey72", colour = "white") +
  geom_vline(data = obs, aes(xintercept = observed_auc),
             colour = "firebrick", linewidth = .9) +
  geom_text(data = obs, aes(x = observed_auc, y = Inf, label = lab),
            hjust = 1.08, vjust = 1.4, size = 2.9, colour = "firebrick") +
  facet_wrap(~phenotype) +
  coord_cartesian(xlim = c(0.2, 1)) +
  labs(x = "ROC-AUC under permuted labels (ridge, 200 permutations)",
       y = "Frequency")
ggsave(file.path(FIG, "fig4_permutation.pdf"), p4, width = 9, height = 3.2)
ggsave(file.path(FIG, "fig4_permutation.png"), p4, width = 9, height = 3.2, dpi = 160)

# ---- Figure 5: ridge penalty plateau --------------------------------------
cv <- read.csv(file.path(RES, "ridge_lambda_curve.csv"))
cv$phenotype <- factor(ph_lab[cv$phenotype], levels = ph_lab)
p5 <- ggplot(cv, aes(lambda, auc, colour = phenotype)) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey55") +
  geom_line(linewidth = .8) + geom_point(size = 1.8) +
  scale_x_log10(breaks = 10^seq(-1, 8, 1),
                labels = c("0.1","1","10","10²","10³","10⁴","10⁵","10⁶","10⁷","10⁸")) +
  scale_colour_brewer(palette = "Dark2") +
  coord_cartesian(ylim = c(0.45, 1)) +
  labs(x = expression("Ridge penalty "*lambda*" (log scale)"),
       y = "ROC-AUC (outer CV)", colour = NULL) +
  theme(legend.position = "bottom")
ggsave(file.path(FIG, "fig5_ridge_plateau.pdf"), p5, width = 6.5, height = 4)
ggsave(file.path(FIG, "fig5_ridge_plateau.png"), p5, width = 6.5, height = 4, dpi = 160)

cat("Figures 4-5 written to", FIG, "\n")
