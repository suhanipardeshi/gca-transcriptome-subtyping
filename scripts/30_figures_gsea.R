# ==============================================================================
# Script: 30_figures_gsea.R
# Purpose: Figures for the GSEA results (script 29).
#   Figure 6 - leading KEGG pathways per phenotype, by normalised enrichment
#              score, coloured by direction.
#   Figure 7 - NES heatmap of pathways significant in ALL THREE phenotypes,
#              which is where the shared immune-up / smooth-muscle-down axis
#              is visible directly.
#   Figure 8 - leading GO biological process terms per phenotype.
# ==============================================================================
setwd("C:/Users/suhan/Desktop/Dissertation_aug")
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
FIG <- "figures/final"; RES <- "results/benchmark"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

ph_lab <- c(intima = "Intima pattern", occlusion = "Occlusion grade",
            media = "Media destruction")

theme_set(theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom"))
# haematoxylin / eosin palette, consistent with the rest of the project
COL <- c(Up = "#C4576B", Down = "#3B2E6E")

read_gsea <- function(kind) {
  do.call(rbind, lapply(names(ph_lab), function(p) {
    f <- file.path(RES, sprintf("gsea_%s_%s.csv", kind, p))
    if (!file.exists(f)) return(NULL)
    d <- read.csv(f, stringsAsFactors = FALSE)
    d$phenotype <- ph_lab[p]
    d[, c("Description", "NES", "p.adjust", "setSize", "phenotype")]
  }))
}
kg <- read_gsea("kegg"); go <- read_gsea("go")
kg$phenotype <- factor(kg$phenotype, levels = ph_lab)
go$phenotype <- factor(go$phenotype, levels = ph_lab)

# ------------------------------------------------------- Figure 6: KEGG bars
top_by_dir <- function(d, n = 8) {
  d %>% group_by(phenotype) %>%
    mutate(dir = ifelse(NES > 0, "Up", "Down")) %>%
    group_by(phenotype, dir) %>%
    slice_max(abs(NES), n = n, with_ties = FALSE) %>% ungroup()
}
k6 <- top_by_dir(kg, 8) %>%
  mutate(lab = ifelse(nchar(Description) > 42,
                      paste0(substr(Description, 1, 40), "\u2026"), Description))

p6 <- ggplot(k6, aes(x = NES, y = reorder(paste(lab, phenotype), NES),
                     fill = dir)) +
  geom_col(width = .72) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = .3) +
  facet_wrap(~phenotype, scales = "free_y", ncol = 3) +
  scale_y_discrete(labels = function(x) sub(" (Intima pattern|Occlusion grade|Media destruction)$", "", x)) +
  scale_fill_manual(values = COL, name = NULL,
                    labels = c(Up = "Up with severity", Down = "Down with severity")) +
  labs(x = "Normalised enrichment score", y = NULL) +
  theme(axis.text.y = element_text(size = 7))
ggsave(file.path(FIG, "fig6_gsea_kegg.pdf"), p6, width = 12, height = 5.6)
ggsave(file.path(FIG, "fig6_gsea_kegg.png"), p6, width = 12, height = 5.6, dpi = 160)

# --------------------------------------- Figure 7: shared-pathway NES heatmap
shared <- kg %>% count(Description) %>% filter(n == 3) %>% pull(Description)
hm <- kg %>% filter(Description %in% shared) %>%
  group_by(Description) %>% mutate(m = mean(abs(NES))) %>% ungroup() %>%
  slice_max(m, n = 3 * 24, with_ties = FALSE) %>%
  mutate(lab = ifelse(nchar(Description) > 46,
                      paste0(substr(Description, 1, 44), "\u2026"), Description))

p7 <- ggplot(hm, aes(phenotype, reorder(lab, m), fill = NES)) +
  geom_tile(colour = "white", linewidth = .4) +
  geom_text(aes(label = sprintf("%+.1f", NES)), size = 2.5,
            colour = ifelse(abs(hm$NES) > 1.7, "white", "grey15")) +
  scale_fill_gradient2(low = COL["Down"], mid = "grey94", high = COL["Up"],
                       midpoint = 0, name = "NES") +
  labs(x = NULL, y = NULL,
       subtitle = sprintf("KEGG pathways significant in all three phenotypes (n = %d shown of %d)",
                          length(unique(hm$Description)), length(shared))) +
  theme(axis.text.y = element_text(size = 7),
        panel.grid = element_blank(),
        plot.subtitle = element_text(size = 9, colour = "grey30"))
ggsave(file.path(FIG, "fig7_gsea_shared.pdf"), p7, width = 8.4, height = 7.4)
ggsave(file.path(FIG, "fig7_gsea_shared.png"), p7, width = 8.4, height = 7.4, dpi = 160)

# ------------------------------------------------------- Figure 8: GO BP dots
g8 <- go %>% group_by(phenotype) %>% slice_max(abs(NES), n = 12, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(dir = ifelse(NES > 0, "Up", "Down"),
         lab = ifelse(nchar(Description) > 44,
                      paste0(substr(Description, 1, 42), "\u2026"), Description))
p8 <- ggplot(g8, aes(NES, reorder(paste(lab, phenotype), NES),
                     colour = dir, size = setSize)) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = .3) +
  geom_point() +
  facet_wrap(~phenotype, scales = "free_y", ncol = 3) +
  scale_y_discrete(labels = function(x) sub(" (Intima pattern|Occlusion grade|Media destruction)$", "", x)) +
  scale_colour_manual(values = COL, name = NULL,
                      labels = c(Up = "Up with severity", Down = "Down with severity")) +
  scale_size_continuous(range = c(1.4, 5), name = "Genes in set") +
  labs(x = "Normalised enrichment score", y = NULL) +
  theme(axis.text.y = element_text(size = 7))
ggsave(file.path(FIG, "fig8_gsea_go.pdf"), p8, width = 12, height = 5.2)
ggsave(file.path(FIG, "fig8_gsea_go.png"), p8, width = 12, height = 5.2, dpi = 160)

cat("KEGG pathways significant in all three phenotypes:", length(shared), "\n")
cat("Figures 6-8 written to", FIG, "\n")
