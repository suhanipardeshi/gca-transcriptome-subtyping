# ==============================================================================
# Script: 24_stability_enrichment.R
# Purpose: Derive gene signatures by stability selection across outer CV folds,
#          quantify signature instability, and run GO over-representation
#          analysis against the correct background universe.
#
# Rationale: a single LASSO fit on one data split yields one gene list and no
# indication of how reproducible that list is. Refitting across the 15 outer
# folds and recording selection frequency measures reproducibility directly.
# ==============================================================================

setwd("C:/Users/suhan/Desktop/Dissertation_aug")
suppressPackageStartupMessages({ library(glmnet) })
RES <- "results/benchmark"; dir.create(RES, showWarnings = FALSE, recursive = TRUE)
N_OUTER <- 5; N_REPEAT <- 3; N_INNER <- 5; SEED <- 42

counts <- as.matrix(readRDS("processed/vst_counts_ml.rds"))
storage.mode(counts) <- "double"
meta <- readRDS("processed/meta_aligned.rds")
num01 <- function(x) as.numeric(as.character(x))

auc_fast <- function(p, y) { n1 <- sum(y==1); n0 <- sum(y==0)
  if (n1==0||n0==0) return(NA_real_); r <- rank(p)
  (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0) }
strat_folds <- function(y, k, seed) { set.seed(seed); i <- rep(NA_integer_, length(y))
  for (cl in c(0,1)) { w <- which(y==cl); i[w] <- sample(rep_len(seq_len(k), length(w))) }; i }
smote_matrix <- function(X, y, k=5, seed=1) { set.seed(seed)
  n1 <- sum(y==1); n0 <- sum(y==0); minc <- if (n1<n0) 1 else 0; need <- abs(n1-n0)
  if (need==0) return(list(X=X,y=y))
  Xi <- X[y==minc,,drop=FALSE]; m <- nrow(Xi); if (m<2) return(list(X=X,y=y))
  kk <- min(k,m-1); D <- as.matrix(dist(Xi)); a <- sample(seq_len(m), need, replace=TRUE)
  Xn <- matrix(0, need, ncol(X))
  for (i in seq_len(need)) { nb <- order(D[a[i],])[2:(kk+1)]; b <- nb[sample.int(kk,1)]
    Xn[i,] <- Xi[a[i],] + runif(1)*(Xi[b,]-Xi[a[i],]) }
  list(X=rbind(X,Xn), y=c(y,rep(minc,need))) }

LAMBDA <- 10^seq(0.5, -4, length.out = 40)
tune_lasso <- function(X, y, seed) {
  fi <- strat_folds(y, N_INNER, seed); sc <- matrix(NA_real_, N_INNER, length(LAMBDA))
  for (f in seq_len(N_INNER)) {
    tr <- fi!=f; te <- fi==f; if (length(unique(y[te]))<2) next
    sm <- smote_matrix(X[tr,,drop=FALSE], y[tr], seed=seed+f)
    m <- glmnet(sm$X, factor(sm$y), family="binomial", alpha=1, lambda=LAMBDA)
    P <- predict(m, newx=X[te,,drop=FALSE], type="response")
    for (j in seq_len(ncol(P))) sc[f,j] <- auc_fast(P[,j], y[te])
  }
  LAMBDA[which.max(colMeans(sc, na.rm=TRUE))] }

PH <- list(
  intima    = list(lab="Intimal inflammation severity", y=num01(meta$Intima_pattern_cat),
                   cc=!is.na(meta$Intima_pattern)),
  occlusion = list(lab="Arterial occlusion grade", y=num01(meta$Occlusion_grade_cat),
                   cc=!is.na(meta$Occlusion_grade)),
  media     = list(lab="Media destruction", y=num01(meta$Media_destruction),
                   cc=!is.na(num01(meta$Media_destruction))))

summ <- list()
for (nm in names(PH)) {
  ph <- PH[[nm]]; keep <- ph$cc & !is.na(ph$y)
  X <- counts[keep,,drop=FALSE]; y <- as.integer(ph$y[keep])
  genes <- colnames(X)
  cat("\n=====", ph$lab, "(n =", length(y), ") =====\n")

  sel <- matrix(FALSE, N_OUTER*N_REPEAT, ncol(X), dimnames=list(NULL, genes))
  for (rep in seq_len(N_REPEAT)) {
    fo <- strat_folds(y, N_OUTER, SEED + 1000*rep)
    for (f in seq_len(N_OUTER)) {
      k <- (rep-1)*N_OUTER + f; tr <- fo!=f
      lam <- tune_lasso(X[tr,,drop=FALSE], y[tr], SEED+k)
      sm  <- smote_matrix(X[tr,,drop=FALSE], y[tr], seed=SEED+k)
      m   <- glmnet(sm$X, factor(sm$y), family="binomial", alpha=1, lambda=LAMBDA)
      co  <- as.matrix(coef(m, s=lam))[-1,1]
      sel[k, names(co)[co != 0]] <- TRUE
    }
  }
  freq <- colSums(sel)
  freq <- freq[freq > 0]
  freq <- sort(freq, decreasing = TRUE)

  cat("genes selected at least once :", length(freq), "\n")
  cat("selected in >=50% of folds   :", sum(freq >= 8), "\n")
  cat("selected in >=80% of folds   :", sum(freq >= 12), "\n")
  cat("selected in ALL 15 folds     :", sum(freq == 15), "\n")
  cat("median folds per selected gene:", median(freq), "\n")

  out <- data.frame(gene = names(freq), n_folds = as.integer(freq),
                    frequency = round(as.numeric(freq)/15, 3))
  write.csv(out, file.path(RES, sprintf("stability_%s.csv", nm)), row.names = FALSE)

  summ[[nm]] <- data.frame(phenotype = nm, label = ph$lab, n = length(y),
    n_ever = length(freq), n_50pc = sum(freq >= 8), n_80pc = sum(freq >= 12),
    n_all = sum(freq == 15), median_folds = median(freq))
}
st <- do.call(rbind, summ)
write.csv(st, file.path(RES, "stability_summary.csv"), row.names = FALSE)
cat("\n===== STABILITY SUMMARY =====\n"); print(st, row.names = FALSE)

# ----------------------------------------------------- GO over-representation
ok <- requireNamespace("clusterProfiler", quietly=TRUE) &&
      requireNamespace("org.Hs.eg.db", quietly=TRUE)
if (ok) {
  suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db)})
  universe <- colnames(counts)          # the 10,000 genes that entered the models
  for (nm in names(PH)) {
    s <- read.csv(file.path(RES, sprintf("stability_%s.csv", nm)))
    gl <- s$gene[s$n_folds >= 8]
    cat("\n--- GO BP:", nm, "| consensus genes (>=50% folds):", length(gl), "---\n")
    if (length(gl) < 5) { cat("too few consensus genes for enrichment\n"); next }
    go <- tryCatch(enrichGO(gene=gl, universe=universe, OrgDb=org.Hs.eg.db,
                            keyType="ENSEMBL", ont="BP", pAdjustMethod="BH",
                            pvalueCutoff=0.05, qvalueCutoff=0.2, readable=TRUE),
                   error=function(e) NULL)
    if (is.null(go) || nrow(as.data.frame(go))==0) {
      cat("no GO BP terms significant at FDR < 0.05 with the correct universe\n")
      write.csv(data.frame(), file.path(RES, sprintf("go_%s.csv", nm)), row.names=FALSE)
    } else {
      d <- as.data.frame(go)
      write.csv(d, file.path(RES, sprintf("go_%s.csv", nm)), row.names=FALSE)
      print(head(d[, c("Description","GeneRatio","pvalue","p.adjust")], 8))
    }
  }
} else cat("\nclusterProfiler/org.Hs.eg.db unavailable; enrichment skipped\n")
cat("\nDONE\n")
