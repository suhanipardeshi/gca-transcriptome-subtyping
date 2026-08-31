# ==============================================================================
# Script: 27_ridge_limit_permutation.R
# Purpose:
#   (a) Show how ridge performance behaves as the penalty grows, to establish
#       that the boundary selection at lambda = 1e6 reflects convergence to the
#       infinite-penalty limit rather than an unbracketed optimum. As lambda
#       increases the coefficient vector shrinks towards zero but its DIRECTION
#       converges; ROC-AUC depends only on the ranking of the linear predictor,
#       which is invariant to positive rescaling, so AUC must plateau.
#   (b) Run a label-permutation test on the ridge model itself. The earlier
#       permutation test in script 22 used the LASSO pipeline, which is the
#       weakest of the six algorithms, so its p-values do not licence any claim
#       about the best-performing model.
# ==============================================================================
setwd("C:/Users/suhan/Desktop/Dissertation_aug")
suppressPackageStartupMessages({ library(glmnet) })
RES <- "results/benchmark"
N_OUTER <- 5; N_REPEAT <- 3; SEED <- 42; N_PERM <- 200
LAMBDA_FIX <- 1e6

counts <- as.matrix(readRDS("processed/vst_counts_ml.rds"))
storage.mode(counts) <- "double"
meta <- readRDS("processed/meta_aligned.rds")
num01 <- function(x) as.numeric(as.character(x))

auc_fast <- function(p, y) { n1 <- sum(y==1); n0 <- sum(y==0)
  if (n1==0||n0==0) return(NA_real_); r <- rank(p)
  (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0) }
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

# outer-CV AUC for ridge at a FIXED lambda (no inner tuning)
ridge_cv <- function(X, y, lambda, seed_base = SEED) {
  a <- numeric()
  for (rep in seq_len(N_REPEAT)) {
    fo <- strat_folds(y, N_OUTER, seed_base + 1000*rep)
    for (f in seq_len(N_OUTER)) {
      k <- (rep-1)*N_OUTER + f; tr <- fo!=f; te <- fo==f
      sm <- smote_matrix(X[tr,,drop=FALSE], y[tr], seed = seed_base + k)
      m  <- glmnet(sm$X, factor(sm$y), family="binomial", alpha=0, lambda=lambda)
      p  <- as.numeric(predict(m, newx=X[te,,drop=FALSE], type="response"))
      a  <- c(a, auc_fast(p, y[te]))
    }
  }
  mean(a, na.rm = TRUE)
}

PH <- list(
  intima    = list(lab="Intima pattern",   y=num01(meta$Intima_pattern_cat),
                   cc=!is.na(meta$Intima_pattern)),
  occlusion = list(lab="Occlusion grade",  y=num01(meta$Occlusion_grade_cat),
                   cc=!is.na(meta$Occlusion_grade)),
  media     = list(lab="Media destruction",y=num01(meta$Media_destruction),
                   cc=!is.na(num01(meta$Media_destruction))))

LAMS <- 10^c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8)
curve <- list(); res <- list(); pdist <- list()

for (nm in names(PH)) {
  ph <- PH[[nm]]; keep <- ph$cc & !is.na(ph$y)
  X <- counts[keep,,drop=FALSE]; y <- as.integer(ph$y[keep])
  cat("\n=====", ph$lab, "(n =", length(y), ") =====\n")

  cat("(a) ridge AUC as a function of penalty:\n")
  aucs <- sapply(LAMS, function(l) ridge_cv(X, y, l))
  for (i in seq_along(LAMS))
    cat(sprintf("    lambda = %-8g  AUC = %.4f\n", LAMS[i], aucs[i]))
  curve[[nm]] <- data.frame(phenotype = nm, lambda = LAMS, auc = aucs)

  obs <- ridge_cv(X, y, LAMBDA_FIX)
  cat(sprintf("(b) observed ridge AUC at lambda = 1e6 : %.4f\n", obs))

  set.seed(SEED + 11)
  perm_labels <- replicate(N_PERM, sample(y), simplify = FALSE)
  stopifnot("permutations not distinct" =
    length(unique(vapply(perm_labels, paste, character(1), collapse=""))) > N_PERM/2)
  pv <- numeric(N_PERM)
  for (i in seq_len(N_PERM)) {
    pv[i] <- ridge_cv(X, perm_labels[[i]], LAMBDA_FIX)
    if (i %% 50 == 0) cat("    permutation", i, "/", N_PERM, "\n")
  }
  p_emp <- (sum(pv >= obs) + 1) / (N_PERM + 1)
  cat(sprintf("    null mean %.3f (95%% range %.3f-%.3f); EMPIRICAL p = %.4f\n",
              mean(pv), quantile(pv,.025), quantile(pv,.975), p_emp))
  res[[nm]] <- data.frame(phenotype=nm, label=ph$lab, n=length(y),
                          observed_auc=obs, null_mean=mean(pv),
                          null_lo=quantile(pv,.025), null_hi=quantile(pv,.975),
                          perm_p=p_emp)
  pdist[[nm]] <- data.frame(phenotype=nm, perm_auc=pv)
}

write.csv(do.call(rbind, curve), file.path(RES,"ridge_lambda_curve.csv"), row.names=FALSE)
write.csv(do.call(rbind, res),   file.path(RES,"ridge_permutation.csv"),  row.names=FALSE)
write.csv(do.call(rbind, pdist), file.path(RES,"ridge_permutation_distribution.csv"), row.names=FALSE)
cat("\n===== RIDGE PERMUTATION SUMMARY =====\n")
print(do.call(rbind, res), row.names = FALSE, digits = 3)
