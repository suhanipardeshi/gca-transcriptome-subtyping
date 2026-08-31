# ==============================================================================
# Script: 21_benchmark_fast.R
# Purpose: Six-algorithm benchmark, three histological phenotypes, under
#          repeated stratified NESTED cross-validation.
#
# Implementation note: models are called directly on numeric matrices rather
# than through tidymodels recipes. With p = 10,000 the recipe/data.frame layer
# dominates runtime; the statistical procedure is unchanged. SMOTE is
# implemented explicitly (Chawla et al., 2002) so that its behaviour in the
# p >> n setting is fully auditable, and is applied ONLY to inner-training
# folds - never to held-out data.
#
# Author: Suhani Balsing Pardeshi
# ==============================================================================

setwd("C:/Users/suhan/Desktop/Dissertation_aug")
suppressPackageStartupMessages({
  library(glmnet); library(ranger); library(xgboost); library(kernlab)
})

RESULT_DIR <- "results/benchmark"
dir.create(RESULT_DIR, showWarnings = FALSE, recursive = TRUE)

N_OUTER <- 5; N_REPEAT <- 3; N_INNER <- 5; SEED <- 42

counts <- as.matrix(readRDS("processed/vst_counts_ml.rds"))
meta   <- readRDS("processed/meta_aligned.rds")
num01  <- function(x) as.numeric(as.character(x))
storage.mode(counts) <- "double"

# ----------------------------------------------------------------- utilities
auc_fast <- function(prob, y01) {                 # Mann-Whitney U statistic
  n1 <- sum(y01 == 1); n0 <- sum(y01 == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(prob)
  (sum(r[y01 == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

pr_auc_fast <- function(prob, y01) {              # average precision
  o <- order(prob, decreasing = TRUE)
  y <- y01[o]
  tp <- cumsum(y); fp <- cumsum(1 - y)
  prec <- tp / (tp + fp); rec <- tp / sum(y)
  sum(diff(c(0, rec)) * prec)
}

strat_folds <- function(y01, k, seed) {           # stratified fold assignment
  set.seed(seed)
  idx <- rep(NA_integer_, length(y01))
  for (cl in c(0, 1)) {
    w <- which(y01 == cl)
    idx[w] <- sample(rep_len(seq_len(k), length(w)))
  }
  idx
}

smote_matrix <- function(X, y01, k = 5, seed = 1) {
  set.seed(seed)
  n1 <- sum(y01 == 1); n0 <- sum(y01 == 0)
  minc <- if (n1 < n0) 1 else 0
  need <- abs(n1 - n0)
  if (need == 0) return(list(X = X, y = y01))
  Xi <- X[y01 == minc, , drop = FALSE]
  m <- nrow(Xi)
  if (m < 2) return(list(X = X, y = y01))
  kk <- min(k, m - 1)
  D  <- as.matrix(dist(Xi))
  a_idx <- sample(seq_len(m), need, replace = TRUE)
  Xnew <- matrix(0, need, ncol(X))
  for (i in seq_len(need)) {
    a  <- a_idx[i]
    nb <- order(D[a, ])[2:(kk + 1)]
    b  <- nb[sample.int(kk, 1)]
    lam <- runif(1)
    Xnew[i, ] <- Xi[a, ] + lam * (Xi[b, ] - Xi[a, ])
  }
  list(X = rbind(X, Xnew), y = c(y01, rep(minc, need)))
}

# ------------------------------------------------------------------- learners
# Each learner: tune on inner folds, refit on full outer-train, return P(Severe)
GRID <- list(
  lasso = list(lambda = 10^seq(0.5, -4, length.out = 40), alpha = 1),
  # Ridge grid extended to 10^6: at an upper bound of 10^3 the optimum was
  # selected AT the boundary for all three phenotypes, meaning the search had
  # not yet bracketed the optimum. Heavy shrinkage is expected here because the
  # signal is diffuse, so the grid must extend far enough to contain the peak.
  ridge = list(lambda = 10^seq(6, -2, length.out = 50), alpha = 0),
  enet  = list(lambda = 10^seq(0.5, -4, length.out = 40),
               alpha  = c(0.15, 0.35, 0.5, 0.7, 0.85)),
  svm   = list(C = c(0.25, 1, 4), sigma = c(1e-4, 1e-3)),
  rf    = list(mtry = c(50, 300), min_node = c(2, 5)),
  # XGBoost search space kept deliberately compact: each fit costs ~3.5 s at
  # p = 10,000, and with 45 outer folds an 8-point grid is not affordable.
  # Learning rates are raised to compensate for the reduced number of rounds.
  xgb   = list(depth = c(2, 4), eta = c(0.1, 0.3), colsample = 0.05)
)

fit_glmnet <- function(Xtr, ytr, Xte, alpha, lambda, seed) {
  sm <- smote_matrix(Xtr, ytr, seed = seed)
  m <- glmnet(sm$X, factor(sm$y), family = "binomial",
              alpha = alpha, lambda = lambda, standardize = TRUE)
  predict(m, newx = Xte, type = "response")          # matrix: n_test x n_lambda
}

tune_glmnet <- function(X, y01, alphas, lambda, seed) {
  fi <- strat_folds(y01, N_INNER, seed)
  best <- list(auc = -Inf)
  for (a in alphas) {
    scores <- matrix(NA_real_, N_INNER, length(lambda))
    for (f in seq_len(N_INNER)) {
      tr <- fi != f; te <- fi == f
      if (length(unique(y01[te])) < 2) next
      P <- fit_glmnet(X[tr, , drop = FALSE], y01[tr], X[te, , drop = FALSE],
                      a, lambda, seed + f)
      for (j in seq_len(ncol(P))) scores[f, j] <- auc_fast(P[, j], y01[te])
    }
    mm <- colMeans(scores, na.rm = TRUE)
    j  <- which.max(mm)
    if (length(j) && is.finite(mm[j]) && mm[j] > best$auc)
      best <- list(auc = mm[j], alpha = a, lambda = lambda[j])
  }
  best
}

learner <- function(name, Xtr, ytr, Xte, seed) {
  g <- GRID[[name]]
  if (name %in% c("lasso", "ridge", "enet")) {
    b <- tune_glmnet(Xtr, ytr, g$alpha, g$lambda, seed)
    P <- fit_glmnet(Xtr, ytr, Xte, b$alpha, b$lambda, seed)
    return(list(prob = as.numeric(P[, 1]),
                param = sprintf("alpha=%.2f lambda=%.5g", b$alpha, b$lambda)))
  }
  # generic inner-CV tuner for the non-glmnet learners
  cfg <- switch(name,
    svm = expand.grid(C = g$C, sigma = g$sigma),
    rf  = expand.grid(mtry = g$mtry, min_node = g$min_node),
    xgb = expand.grid(depth = g$depth, eta = g$eta, colsample = g$colsample))
  fi <- strat_folds(ytr, N_INNER, seed)
  # SMOTE depends only on the fold and seed, not on the hyperparameter config,
  # so expand each inner-training fold ONCE and reuse it across all configs.
  cache <- vector("list", N_INNER)
  for (f in seq_len(N_INNER)) {
    tr <- fi != f
    cache[[f]] <- smote_matrix(Xtr[tr, , drop = FALSE], ytr[tr], seed = seed + f)
  }
  sc <- rep(NA_real_, nrow(cfg))
  for (ci in seq_len(nrow(cfg))) {
    vals <- rep(NA_real_, N_INNER)
    for (f in seq_len(N_INNER)) {
      te <- fi == f
      if (length(unique(ytr[te])) < 2) next
      vals[f] <- tryCatch(
        auc_fast(raw_fit(name, cache[[f]]$X, cache[[f]]$y,
                         Xtr[te, , drop = FALSE], cfg[ci, ], seed + f), ytr[te]),
        error = function(e) NA_real_)
    }
    sc[ci] <- mean(vals, na.rm = TRUE)
  }
  ci <- which.max(sc)
  if (!length(ci)) ci <- 1
  smf <- smote_matrix(Xtr, ytr, seed = seed)
  list(prob = raw_fit(name, smf$X, smf$y, Xte, cfg[ci, ], seed),
       param = paste(names(cfg), unlist(cfg[ci, ]), sep = "=", collapse = " "))
}

# NOTE: raw_fit now receives ALREADY SMOTE-expanded training data.
raw_fit <- function(name, X, y, Xte, cfg, seed) {
  if (name == "svm") {
    # Pre-scale once and disable kernlab's internal scaling (which would
    # rescale all 10,000 columns on every fit). Platt probability scaling is
    # also disabled: it runs a further internal CV per fit and is unnecessary,
    # because ROC-AUC depends only on the ranking of the decision values.
    mu <- colMeans(X); sdv <- apply(X, 2, sd); sdv[sdv == 0] <- 1
    Xs  <- scale(X,   center = mu, scale = sdv)
    Xts <- scale(Xte, center = mu, scale = sdv)
    m <- kernlab::ksvm(Xs, factor(y), type = "C-svc", kernel = "rbfdot",
                       kpar = list(sigma = cfg$sigma), C = cfg$C,
                       prob.model = FALSE, scaled = FALSE)
    dtr <- as.numeric(kernlab::predict(m, Xs,  type = "decision"))
    dte <- as.numeric(kernlab::predict(m, Xts, type = "decision"))
    # orient the decision value so that larger = more likely "Severe" (y = 1)
    sgn <- if (mean(dtr[y == 1]) < mean(dtr[y == 0])) -1 else 1
    return(stats::plogis(dte * sgn))   # monotone map; preserves AUC, 0 -> 0.5
  }
  if (name == "rf") {
    m <- ranger::ranger(x = X, y = factor(y), num.trees = 500,
                        mtry = cfg$mtry, min.node.size = cfg$min_node,
                        probability = TRUE, num.threads = 6, seed = seed)
    return(predict(m, Xte)$predictions[, "1"])
  }
  if (name == "xgb") {
    # xgboost >= 3.0 changed the xgboost() signature; xgb.train + DMatrix is
    # the stable interface across versions.
    dtr <- xgboost::xgb.DMatrix(data = X, label = as.numeric(y))
    m <- xgboost::xgb.train(
      params = list(max_depth = cfg$depth, eta = cfg$eta,
                    colsample_bytree = cfg$colsample,
                    objective = "binary:logistic", eval_metric = "logloss",
                    tree_method = "hist", nthread = 2),
      data = dtr, nrounds = 100, verbose = 0)
    return(as.numeric(predict(m, xgboost::xgb.DMatrix(data = Xte))))
  }
  stop("unknown learner")
}

# ---------------------------------------------------------------- phenotypes
PH <- list(
  intima    = list(label = "Intimal inflammation severity",
                   y = num01(meta$Intima_pattern_cat),
                   cc = !is.na(meta$Intima_pattern)),
  occlusion = list(label = "Arterial occlusion grade",
                   y = num01(meta$Occlusion_grade_cat),
                   cc = !is.na(meta$Occlusion_grade)),
  media     = list(label = "Media destruction",
                   y = num01(meta$Media_destruction),
                   cc = !is.na(num01(meta$Media_destruction)))
)

run_one <- function(ph_name, algo) {
  ph <- PH[[ph_name]]
  keep <- ph$cc & !is.na(ph$y)
  X <- counts[keep, , drop = FALSE]
  y <- as.integer(ph$y[keep])            # 1 = Severe

  rows <- list(); allp <- list(); prm <- character()
  for (rep in seq_len(N_REPEAT)) {
    fo <- strat_folds(y, N_OUTER, SEED + 1000 * rep)
    for (f in seq_len(N_OUTER)) {
      tr <- fo != f; te <- fo == f
      k <- (rep - 1) * N_OUTER + f
      r <- learner(algo, X[tr, , drop = FALSE], y[tr],
                   X[te, , drop = FALSE], SEED + k)
      pr <- r$prob; yt <- y[te]
      cl <- as.integer(pr >= 0.5)
      tp <- sum(cl == 1 & yt == 1); tn <- sum(cl == 0 & yt == 0)
      fp <- sum(cl == 1 & yt == 0); fn <- sum(cl == 0 & yt == 1)
      sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA
      spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA
      prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA
      f1   <- if (!is.na(prec) && !is.na(sens) && (prec + sens) > 0)
                2 * prec * sens / (prec + sens) else NA
      rows[[k]] <- data.frame(fold = k, roc_auc = auc_fast(pr, yt),
                              pr_auc = pr_auc_fast(pr, yt),
                              bal_accuracy = mean(c(sens, spec), na.rm = TRUE),
                              sensitivity = sens, specificity = spec, f1 = f1)
      allp[[k]] <- data.frame(fold = k, prob = pr, truth = yt)
      prm <- c(prm, r$param)
      cat(sprintf("    fold %2d  AUC=%.3f  [%s]\n", k, rows[[k]]$roc_auc, r$param))
    }
  }
  ft <- do.call(rbind, rows); ap <- do.call(rbind, allp)
  write.csv(ft, file.path(RESULT_DIR, sprintf("folds_%s_%s.csv", ph_name, algo)), row.names = FALSE)
  write.csv(ap, file.path(RESULT_DIR, sprintf("preds_%s_%s.csv", ph_name, algo)), row.names = FALSE)
  s <- data.frame(
    phenotype = ph_name, phenotype_label = ph$label, algorithm = algo,
    n = length(y), n_severe = sum(y == 1), n_mild = sum(y == 0),
    n_folds = nrow(ft),
    roc_auc_mean = mean(ft$roc_auc, na.rm = TRUE),
    roc_auc_sd   = sd(ft$roc_auc, na.rm = TRUE),
    roc_auc_pooled = auc_fast(ap$prob, ap$truth),
    pr_auc_mean = mean(ft$pr_auc, na.rm = TRUE),
    bal_acc_mean = mean(ft$bal_accuracy, na.rm = TRUE),
    sens_mean = mean(ft$sensitivity, na.rm = TRUE),
    spec_mean = mean(ft$specificity, na.rm = TRUE),
    f1_mean = mean(ft$f1, na.rm = TRUE),
    modal_param = names(sort(table(prm), decreasing = TRUE))[1])
  write.csv(s, file.path(RESULT_DIR, sprintf("summary_%s_%s.csv", ph_name, algo)), row.names = FALSE)
  s
}

args <- commandArgs(trailingOnly = TRUE)
phs <- if (length(args) >= 1) strsplit(args[1], ",")[[1]] else names(PH)
als <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else
         c("lasso", "ridge", "enet", "svm", "rf", "xgb")

for (p in phs) for (a in als) {
  fn <- file.path(RESULT_DIR, sprintf("summary_%s_%s.csv", p, a))
  if (file.exists(fn)) { cat(sprintf("\n=== %s | %s === (done, skip)\n", p, a)); next }
  cat(sprintf("\n=== %s | %s ===\n", PH[[p]]$label, a))
  t0 <- Sys.time()
  s <- try(run_one(p, a), silent = TRUE)
  if (inherits(s, "try-error")) cat("  FAILED:", as.character(s), "\n")
  else cat(sprintf("  -> AUC %.3f (SD %.3f) pooled %.3f  [%.1f min]\n",
                   s$roc_auc_mean, s$roc_auc_sd, s$roc_auc_pooled,
                   as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  flush.console()
}
cat("\nDONE\n")
