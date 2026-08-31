# ==============================================================================
# Script: 22_baselines_nulls.R
# Purpose: Establish what the transcriptomic models must beat, and test whether
#          the observed performance exceeds chance.
#
#   A. Simple baselines  - sequencing depth alone; clinical variables alone;
#                          depth + clinical. Evaluated under the SAME repeated
#                          stratified CV as the transcriptomic models.
#   B. Random-signature null - LASSO restricted to k randomly chosen genes,
#                          k = size of the real signature. Tests whether the
#                          selected genes outperform arbitrary genes.
#   C. Label-permutation null - outcome labels shuffled, full nested CV re-run.
#                          Yields an empirical p-value for the observed AUC.
#   D. Depth-adjusted models - expression residualised on log sequencing depth
#                          before modelling, to test whether signal survives
#                          removal of the technical axis.
#
# Author: Suhani Balsing Pardeshi
# ==============================================================================

setwd("C:/Users/suhan/Desktop/Dissertation_aug")
suppressPackageStartupMessages({ library(glmnet) })

RESULT_DIR <- "results/benchmark"
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)
N_OUTER <- 5; N_REPEAT <- 3; N_INNER <- 5; SEED <- 42
N_PERM  <- 200; N_RANDSIG <- 200

counts <- as.matrix(readRDS("processed/vst_counts_ml.rds"))
meta   <- readRDS("processed/meta_aligned.rds")
storage.mode(counts) <- "double"
num01  <- function(x) as.numeric(as.character(x))

auc_fast <- function(prob, y01) {
  n1 <- sum(y01 == 1); n0 <- sum(y01 == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(prob); (sum(r[y01 == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
strat_folds <- function(y01, k, seed) {
  set.seed(seed); idx <- rep(NA_integer_, length(y01))
  for (cl in c(0, 1)) { w <- which(y01 == cl)
    idx[w] <- sample(rep_len(seq_len(k), length(w))) }
  idx
}
smote_matrix <- function(X, y01, k = 5, seed = 1) {
  set.seed(seed)
  n1 <- sum(y01 == 1); n0 <- sum(y01 == 0)
  minc <- if (n1 < n0) 1 else 0; need <- abs(n1 - n0)
  if (need == 0) return(list(X = X, y = y01))
  Xi <- X[y01 == minc, , drop = FALSE]; m <- nrow(Xi)
  if (m < 2) return(list(X = X, y = y01))
  kk <- min(k, m - 1); D <- as.matrix(dist(Xi))
  a_idx <- sample(seq_len(m), need, replace = TRUE)
  Xnew <- matrix(0, need, ncol(X))
  for (i in seq_len(need)) {
    a <- a_idx[i]; nb <- order(D[a, ])[2:(kk + 1)]
    b <- nb[sample.int(kk, 1)]; lam <- runif(1)
    Xnew[i, ] <- Xi[a, ] + lam * (Xi[b, ] - Xi[a, ])
  }
  list(X = rbind(X, Xnew), y = c(y01, rep(minc, need)))
}

LAMBDA <- 10^seq(0.5, -4, length.out = 40)
tune_lasso <- function(X, y01, seed) {
  fi <- strat_folds(y01, N_INNER, seed)
  sc <- matrix(NA_real_, N_INNER, length(LAMBDA))
  for (f in seq_len(N_INNER)) {
    tr <- fi != f; te <- fi == f
    if (length(unique(y01[te])) < 2) next
    sm <- smote_matrix(X[tr, , drop = FALSE], y01[tr], seed = seed + f)
    m <- glmnet(sm$X, factor(sm$y), family = "binomial", alpha = 1, lambda = LAMBDA)
    P <- predict(m, newx = X[te, , drop = FALSE], type = "response")
    for (j in seq_len(ncol(P))) sc[f, j] <- auc_fast(P[, j], y01[te])
  }
  mm <- colMeans(sc, na.rm = TRUE); LAMBDA[which.max(mm)]
}

nested_lasso_auc <- function(X, y, seed_base = SEED) {
  ap <- list(); aucs <- numeric()
  for (rep in seq_len(N_REPEAT)) {
    fo <- strat_folds(y, N_OUTER, seed_base + 1000 * rep)
    for (f in seq_len(N_OUTER)) {
      tr <- fo != f; te <- fo == f
      k <- (rep - 1) * N_OUTER + f
      lam <- tune_lasso(X[tr, , drop = FALSE], y[tr], seed_base + k)
      sm <- smote_matrix(X[tr, , drop = FALSE], y[tr], seed = seed_base + k)
      m <- glmnet(sm$X, factor(sm$y), family = "binomial", alpha = 1, lambda = LAMBDA)
      p <- as.numeric(predict(m, newx = X[te, , drop = FALSE],
                              type = "response", s = lam))
      aucs <- c(aucs, auc_fast(p, y[te]))
      ap[[k]] <- data.frame(prob = p, truth = y[te])
    }
  }
  all <- do.call(rbind, ap)
  list(mean = mean(aucs, na.rm = TRUE), sd = sd(aucs, na.rm = TRUE),
       pooled = auc_fast(all$prob, all$truth))
}

# simple logistic baseline under identical CV structure
nested_glm_auc <- function(Z, y, seed_base = SEED) {
  ap <- list(); aucs <- numeric()
  Z <- as.data.frame(Z)
  for (rep in seq_len(N_REPEAT)) {
    fo <- strat_folds(y, N_OUTER, seed_base + 1000 * rep)
    for (f in seq_len(N_OUTER)) {
      tr <- fo != f; te <- fo == f
      k <- (rep - 1) * N_OUTER + f
      d <- data.frame(y = y[tr], Z[tr, , drop = FALSE])
      fit <- suppressWarnings(glm(y ~ ., data = d, family = binomial))
      p <- suppressWarnings(predict(fit, newdata = Z[te, , drop = FALSE],
                                    type = "response"))
      aucs <- c(aucs, auc_fast(p, y[te]))
      ap[[k]] <- data.frame(prob = p, truth = y[te])
    }
  }
  all <- do.call(rbind, ap)
  list(mean = mean(aucs, na.rm = TRUE), sd = sd(aucs, na.rm = TRUE),
       pooled = auc_fast(all$prob, all$truth))
}

PH <- list(
  intima    = list(label = "Intimal inflammation severity",
                   y = num01(meta$Intima_pattern_cat), cc = !is.na(meta$Intima_pattern),
                   sig_k = 38),
  occlusion = list(label = "Arterial occlusion grade",
                   y = num01(meta$Occlusion_grade_cat), cc = !is.na(meta$Occlusion_grade),
                   sig_k = 33),
  media     = list(label = "Media destruction",
                   y = num01(meta$Media_destruction),
                   cc = !is.na(num01(meta$Media_destruction)), sig_k = 34)
)

out <- list(); pn <- list()

for (nm in names(PH)) {
  ph <- PH[[nm]]
  keep <- ph$cc & !is.na(ph$y)
  X <- counts[keep, , drop = FALSE]
  y <- as.integer(ph$y[keep])
  md <- meta[keep, ]
  depth <- log10(md$sequencing_depth)

  cat("\n=====================================================\n")
  cat(ph$label, " (n =", length(y), "; severe =", sum(y), ")\n")
  cat("=====================================================\n")

  # ---- A. simple baselines -------------------------------------------------
  b_depth <- nested_glm_auc(data.frame(depth = depth), y)
  cat(sprintf("A1 depth only          AUC = %.3f (SD %.3f) pooled %.3f\n",
              b_depth$mean, b_depth$sd, b_depth$pooled))

  clin <- data.frame(age = md$age, sex = md$sex,
                     steroids = ifelse(is.na(md$steroids_days),
                                       median(md$steroids_days, na.rm = TRUE),
                                       md$steroids_days))
  b_clin <- nested_glm_auc(clin, y)
  cat(sprintf("A2 clinical only       AUC = %.3f (SD %.3f) pooled %.3f\n",
              b_clin$mean, b_clin$sd, b_clin$pooled))

  b_both <- nested_glm_auc(cbind(clin, depth = depth), y)
  cat(sprintf("A3 clinical + depth    AUC = %.3f (SD %.3f) pooled %.3f\n",
              b_both$mean, b_both$sd, b_both$pooled))

  # ---- D. depth-adjusted expression ---------------------------------------
  # Residualise every gene on log10 depth. This is ordinary least squares with
  # a single predictor, so it is computed in closed form rather than by 10,000
  # separate lm() calls (identical result, orders of magnitude faster).
  dc <- depth - mean(depth)
  Xc <- sweep(X, 2, colMeans(X), "-")
  b  <- as.numeric(crossprod(dc, Xc)) / sum(dc^2)
  Xr <- Xc - outer(dc, b)
  real <- nested_lasso_auc(X, y)
  adj  <- nested_lasso_auc(Xr, y)
  cat(sprintf("D1 LASSO (raw VST)     AUC = %.3f (SD %.3f) pooled %.3f\n",
              real$mean, real$sd, real$pooled))
  cat(sprintf("D2 LASSO depth-adjust  AUC = %.3f (SD %.3f) pooled %.3f\n",
              adj$mean, adj$sd, adj$pooled))

  # ---- B. random-signature null -------------------------------------------
  # NOTE: nested_lasso_auc() calls set.seed() internally (for fold assignment
  # and SMOTE), which resets the global RNG. Random draws must therefore be
  # generated UP FRONT, or every iteration would reuse an identical draw.
  set.seed(SEED)
  rand_cols <- replicate(N_RANDSIG, sample.int(ncol(X), ph$sig_k),
                         simplify = FALSE)
  rs <- numeric(N_RANDSIG)
  for (i in seq_len(N_RANDSIG)) {
    rs[i] <- nested_lasso_auc(X[, rand_cols[[i]], drop = FALSE], y,
                              seed_base = SEED)$mean
  }
  cat(sprintf("B  random %d-gene sig  AUC = %.3f (95%% range %.3f-%.3f)\n",
              ph$sig_k, mean(rs), quantile(rs, .025), quantile(rs, .975)))

  # ---- C. label-permutation null ------------------------------------------
  # Permutations are likewise pre-generated (see note above) so that each
  # iteration uses a genuinely different shuffle of the outcome labels.
  set.seed(SEED + 7)
  perm_labels <- replicate(N_PERM, sample(y), simplify = FALSE)
  stopifnot("permutations are not distinct - RNG was reset" =
              length(unique(vapply(perm_labels, paste, character(1),
                                   collapse = ""))) > N_PERM / 2)
  perm <- numeric(N_PERM)
  for (i in seq_len(N_PERM)) {
    perm[i] <- nested_lasso_auc(X, perm_labels[[i]], seed_base = SEED)$mean
    if (i %% 25 == 0) cat("   permutation", i, "/", N_PERM, "\n")
  }
  p_emp <- (sum(perm >= real$mean) + 1) / (N_PERM + 1)
  cat(sprintf("C  permutation null    AUC = %.3f (95%% range %.3f-%.3f)\n",
              mean(perm), quantile(perm, .025), quantile(perm, .975)))
  cat(sprintf("   EMPIRICAL p-value for observed AUC %.3f : p = %.4f\n",
              real$mean, p_emp))

  out[[nm]] <- data.frame(
    phenotype = nm, label = ph$label, n = length(y), n_severe = sum(y),
    depth_auc = b_depth$mean, clin_auc = b_clin$mean, clin_depth_auc = b_both$mean,
    lasso_auc = real$mean, lasso_sd = real$sd,
    lasso_depth_adj_auc = adj$mean,
    rand_sig_auc = mean(rs), rand_sig_lo = quantile(rs, .025),
    rand_sig_hi = quantile(rs, .975),
    perm_auc = mean(perm), perm_hi95 = quantile(perm, .95), perm_p = p_emp)
  pn[[nm]] <- data.frame(phenotype = nm, perm_auc = perm)
}

res <- do.call(rbind, out)
write.csv(res, file.path(RESULT_DIR, "baselines_and_nulls.csv"), row.names = FALSE)
write.csv(do.call(rbind, pn), file.path(RESULT_DIR, "permutation_distribution.csv"),
          row.names = FALSE)
cat("\n\n================ CONSOLIDATED ================\n")
print(res, row.names = FALSE, digits = 3)
cat("\nSaved to", RESULT_DIR, "\n")
