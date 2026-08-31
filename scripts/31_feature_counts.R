# ==============================================================================
# Script: 31_feature_counts.R
# Purpose: How many genes does each of the six models actually use?
#
# "Number of genes selected" is only well defined for the sparse estimators.
# For the others the comparable quantity is the number of genes that carry any
# influence on the fitted model:
#   LASSO / Elastic net / Ridge - coefficients not exactly zero at the chosen
#                                 penalty (ridge is dense by construction).
#   Random forest               - genes used in at least one split, taken as
#                                 impurity importance > 0.
#   XGBoost                     - genes appearing in the fitted booster, from
#                                 xgb.importance().
#   SVM (RBF)                   - the kernel is computed over all predictors, so
#                                 every gene enters; there is no selection.
#
# Models are refit on each of the 15 outer training folds using the modal
# hyperparameters from the benchmark (script 21), so the counts describe the
# same configurations that produced the reported performance.
#
# Note on the LASSO count specifically: this script holds one fixed lambda per
# phenotype (the modal value from script 21) across all 15 refits. Script
# 24_stability_enrichment.R instead re-tunes lambda within each outer fold's
# own inner cross-validation loop, matching the nested cross-validation design
# in Chapter 3 of the dissertation. The two scripts therefore report different
# LASSO gene-per-fold counts by design. The dissertation (Table 3, Section 4.4)
# quotes script 24's figures, not this script's, because 24 uses the same
# fold-by-fold tuning as the rest of the benchmark.
# ==============================================================================

setwd("C:/Users/suhan/Desktop/Dissertation_aug")
suppressPackageStartupMessages({library(glmnet); library(ranger); library(xgboost)})
RES <- "results/benchmark"
N_OUTER <- 5; N_REPEAT <- 3; SEED <- 42

counts <- as.matrix(readRDS("processed/vst_counts_ml.rds"))
storage.mode(counts) <- "double"
meta <- readRDS("processed/meta_aligned.rds")
num01 <- function(x) as.numeric(as.character(x))
P <- ncol(counts)

strat_folds <- function(y,k,seed){ set.seed(seed); i <- rep(NA_integer_,length(y))
  for(cl in c(0,1)){w<-which(y==cl); i[w]<-sample(rep_len(seq_len(k),length(w)))}; i }
smote_matrix <- function(X,y,k=5,seed=1){ set.seed(seed)
  n1<-sum(y==1);n0<-sum(y==0);minc<-if(n1<n0)1 else 0;need<-abs(n1-n0)
  if(need==0) return(list(X=X,y=y))
  Xi<-X[y==minc,,drop=FALSE];m<-nrow(Xi); if(m<2) return(list(X=X,y=y))
  kk<-min(k,m-1);D<-as.matrix(dist(Xi));a<-sample(seq_len(m),need,replace=TRUE)
  Xn<-matrix(0,need,ncol(X))
  for(i in seq_len(need)){nb<-order(D[a[i],])[2:(kk+1)];b<-nb[sample.int(kk,1)]
    Xn[i,]<-Xi[a[i],]+runif(1)*(Xi[b,]-Xi[a[i],])}
  list(X=rbind(X,Xn),y=c(y,rep(minc,need))) }

# modal hyperparameters from the benchmark
MP <- list(
  intima    = list(lasso=list(a=1,l=0.2219), enet=list(a=0.15,l=1.4251),
                   ridge=list(a=0,l=1e6), rf=list(mtry=50,mn=2),
                   xgb=list(d=4,eta=0.1,cs=0.05)),
  occlusion = list(lasso=list(a=1,l=0.17013), enet=list(a=0.15,l=1.0926),
                   ridge=list(a=0,l=1e6), rf=list(mtry=50,mn=2),
                   xgb=list(d=2,eta=0.1,cs=0.05)),
  media     = list(lasso=list(a=1,l=0.1), enet=list(a=0.15,l=1.4251),
                   ridge=list(a=0,l=1e6), rf=list(mtry=50,mn=2),
                   xgb=list(d=2,eta=0.3,cs=0.05)))

PH <- list(
  intima    = list(lab="Intima pattern",   y=num01(meta$Intima_pattern_cat),
                   cc=!is.na(meta$Intima_pattern)),
  occlusion = list(lab="Occlusion grade",  y=num01(meta$Occlusion_grade_cat),
                   cc=!is.na(meta$Occlusion_grade)),
  media     = list(lab="Media destruction",y=num01(meta$Media_destruction),
                   cc=!is.na(num01(meta$Media_destruction))))

rows <- list(); unions <- list()
for (nm in names(PH)) {
  ph <- PH[[nm]]; keep <- ph$cc & !is.na(ph$y)
  X <- counts[keep,,drop=FALSE]; y <- as.integer(ph$y[keep])
  mp <- MP[[nm]]
  per <- list(lasso=integer(), enet=integer(), ridge=integer(),
              svm=integer(), rf=integer(), xgb=integer())
  uni <- list(lasso=character(), enet=character(), ridge=character(),
              svm=character(), rf=character(), xgb=character())

  for (rep in seq_len(N_REPEAT)) {
    fo <- strat_folds(y, N_OUTER, SEED + 1000*rep)
    for (f in seq_len(N_OUTER)) {
      k <- (rep-1)*N_OUTER + f; tr <- fo != f
      sm <- smote_matrix(X[tr,,drop=FALSE], y[tr], seed = SEED + k)

      for (a in c("lasso","enet","ridge")) {
        g <- glmnet(sm$X, factor(sm$y), family="binomial",
                    alpha = mp[[a]]$a, lambda = mp[[a]]$l)
        co <- as.matrix(coef(g))[-1,1]
        nz <- names(co)[co != 0]
        per[[a]] <- c(per[[a]], length(nz)); uni[[a]] <- union(uni[[a]], nz)
      }

      rf <- ranger(x = sm$X, y = factor(sm$y), num.trees = 500,
                   mtry = mp$rf$mtry, min.node.size = mp$rf$mn,
                   probability = TRUE, importance = "impurity",
                   num.threads = 6, seed = SEED + k)
      used <- names(rf$variable.importance)[rf$variable.importance > 0]
      per$rf <- c(per$rf, length(used)); uni$rf <- union(uni$rf, used)

      dtr <- xgb.DMatrix(data = sm$X, label = as.numeric(sm$y))
      bst <- xgb.train(params = list(max_depth=mp$xgb$d, eta=mp$xgb$eta,
                        colsample_bytree=mp$xgb$cs, objective="binary:logistic",
                        eval_metric="logloss", tree_method="hist", nthread=2),
                       data = dtr, nrounds = 100, verbose = 0)
      imp <- tryCatch(xgb.importance(model = bst), error=function(e) NULL)
      fx <- if (is.null(imp)) character() else imp$Feature
      per$xgb <- c(per$xgb, length(fx)); uni$xgb <- union(uni$xgb, fx)

      per$svm <- c(per$svm, P)   # RBF kernel uses every predictor
    }
    cat(".")
  }
  uni$svm <- colnames(X)
  cat(" ", ph$lab, "\n")

  for (a in c("lasso","enet","ridge","svm","rf","xgb")) {
    v <- per[[a]]
    rows[[length(rows)+1]] <- data.frame(
      phenotype = nm, label = ph$lab, algorithm = a,
      mean_per_fold = mean(v), min_per_fold = min(v), max_per_fold = max(v),
      union_across_folds = length(uni[[a]]),
      pct_of_10000 = round(100*mean(v)/P, 1))
  }
  unions[[nm]] <- uni
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path(RES, "feature_counts.csv"), row.names = FALSE)

cat("\n=========== GENES ACTUALLY USED BY EACH MODEL ===========\n")
alab <- c(lasso="LASSO", enet="Elastic net", ridge="Ridge",
          svm="SVM (RBF)", rf="Random forest", xgb="XGBoost")
for (nm in names(PH)) {
  cat("\n---", PH[[nm]]$lab, "---\n")
  cat(sprintf("%-14s %12s %14s %14s\n", "model", "per fold", "range", "union/15 folds"))
  d <- tab[tab$phenotype == nm, ]
  for (i in seq_len(nrow(d))) {
    cat(sprintf("%-14s %12.1f %7d-%-6d %14d\n", alab[d$algorithm[i]],
                d$mean_per_fold[i], d$min_per_fold[i], d$max_per_fold[i],
                d$union_across_folds[i]))
  }
}
cat("\nSaved to", file.path(RES, "feature_counts.csv"), "\n")
