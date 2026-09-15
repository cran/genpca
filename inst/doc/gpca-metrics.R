params <-
list(family = "red", preset = "interaction")

## ----setup, include = FALSE---------------------------------------------------
if (requireNamespace("ragg", quietly = TRUE)) knitr::opts_chunk$set(dev = "ragg_png")
if (requireNamespace("systemfonts", quietly = TRUE) && requireNamespace("albersdown", quietly = TRUE)) albersdown::albers_register_fonts()
if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("albersdown", quietly = TRUE)) ggplot2::theme_set(albersdown::theme_albers(family = params$family, preset = params$preset))
knitr::opts_chunk$set(
  collapse   = TRUE,
  comment    = "#>",
  message    = FALSE,
  warning    = TRUE,
  fig.width  = 6,
  fig.height = 4,
  out.width  = "85%"
)
library(genpca)
library(Matrix)

flip <- function(M) as.matrix(M)[, ncol(M):1]
heat <- function(M, main = "") {
  image(flip(M), axes = FALSE, main = main,
        col = grey.colors(20, start = 0.95, end = 0.2))
}

## ----albers-classes, echo=FALSE, results='asis'-------------------------------
cat(sprintf(
  paste0(
    '<script>document.addEventListener("DOMContentLoaded",function(){',
    'document.body.classList.remove("palette-red","palette-lapis","palette-ochre","palette-teal","palette-green","palette-violet","preset-homage","preset-interaction","preset-study","preset-structural","preset-adobe","preset-midnight");',
    'document.body.classList.add("palette-%s","preset-%s");',
    '});</script>'
  ),
  params$family,
  params$preset
))

## ----hetero-demo--------------------------------------------------------------
set.seed(42)
n <- 60; p <- 20
X <- matrix(rnorm(n * p), n, p)

col_noise_sd <- runif(p, 0.5, 2)
A <- Diagonal(x = 1 / col_noise_sd^2)
row_noise_sd <- runif(n, 0.7, 1.3)
M <- Diagonal(x = 1 / row_noise_sd^2)

fit <- genpca(X, M = M, A = A, ncomp = 3,
              preproc = multivarious::center())
fit$sdev

## ----hetero-plot, echo = FALSE, fig.cap = "Inverse-variance weights on columns (top) and rows (bottom). Noisier dimensions get smaller weights.", fig.height = 3.5----
op <- par(mfrow = c(2, 1), mar = c(2.5, 4, 2, 1))
barplot(diag(A), border = NA, col = "steelblue",
        main = "Column weights diag(A)", names.arg = NA)
barplot(diag(M), border = NA, col = "tomato",
        main = "Row weights diag(M)", names.arg = NA)
par(op)

## ----scaled-pca-equiv---------------------------------------------------------
set.seed(42)
Xv <- matrix(rnorm(60 * 20), 60, 20) %*% diag(runif(20, 0.5, 3))
sds <- apply(Xv, 2, sd)

g  <- genpca(Xv, A = Diagonal(x = 1 / sds^2), ncomp = 5,
             preproc = multivarious::center())
pr <- prcomp(Xv, scale. = TRUE)

# scores agree component by component
sapply(1:3, function(k) cor(multivarious::scores(g)[, k], pr$x[, k]))

## ----scaled-pca-sdev----------------------------------------------------------
rbind(genpca = g$sdev[1:5],
      prcomp = pr$sdev[1:5],
      ratio  = g$sdev[1:5] / pr$sdev[1:5])
sqrt(nrow(Xv) - 1)

## ----both-margins-------------------------------------------------------------
Xc <- scale(Xv, center = TRUE, scale = FALSE)
W  <- diag(1 / apply(Xc, 1, sd)) %*% Xc %*% diag(1 / apply(Xc, 2, sd))

range(apply(W, 1, sd))   # row SDs, would be constant if standardised
range(apply(W, 2, sd))   # column SDs

## ----scale-indeterminacy------------------------------------------------------
M0 <- Diagonal(x = 1 / apply(Xv, 1, sd)^2)
A0 <- Diagonal(x = 1 / sds^2)

f1 <- genpca(Xv, M = M0,     A = A0,     ncomp = 4, preproc = multivarious::center())
f2 <- genpca(Xv, M = 7 * M0, A = A0 / 7, ncomp = 4, preproc = multivarious::center())

max(abs(f1$sdev - f2$sdev))

## ----loadings-are-AV----------------------------------------------------------
set.seed(1)
Xd <- matrix(rnorm(400), 40, 10)
Ad <- crossprod(matrix(rnorm(100), 10, 10)) / 10 + diag(10)
fd <- genpca(Xd, A = Ad, ncomp = 3, preproc = multivarious::center())

max(abs(multivarious::components(fd) - as.matrix(Ad %*% fd$ov)))

## ----reciprocal-weights-------------------------------------------------------
p <- 32
Wc <- matrix(0, p, p)
for (i in 1:p) { Wc[i, i %% p + 1] <- 1; Wc[i %% p + 1, i] <- 1 }
Lc <- diag(rowSums(Wc)) - Wc

ec <- eigen(Lc, symmetric = TRUE)
o  <- order(ec$values)
Vc <- ec$vectors[, o]          # smoothest first
lc <- ec$values[o]             # Dirichlet energy = roughness

Ac  <- diag(p) + 0.45 * Wc     # a smoother (PSD)
Sig <- solve(Ac)               # the noise covariance it implies

modes <- c(1, 16, 32)          # smoothest, middling, roughest
data.frame(
  roughness   = round(lc[modes], 3),
  metric_wt   = round(sapply(modes, function(k) t(Vc[, k]) %*% Ac  %*% Vc[, k]), 3),
  noise_var   = round(sapply(modes, function(k) t(Vc[, k]) %*% Sig %*% Vc[, k]), 3)
)

## ----orientation-demo, echo = FALSE, fig.cap = "The same graph, two metrics. Left: a smoother concentrates PC1 on the smooth blob. Right: the Laplacian whitens the smooth field away, so PC1 locks onto the fine-scale checkerboard that was buried underneath it.", fig.width = 7, fig.height = 3.5----
set.seed(7)
gg <- 12; pp <- gg * gg; nn <- 80
gidx <- expand.grid(r = 1:gg, c = 1:gg)
Wg <- matrix(0, pp, pp)
for (i in 1:pp) for (j in 1:pp)
  if (i < j && abs(gidx$r[i] - gidx$r[j]) + abs(gidx$c[i] - gidx$c[j]) == 1) {
    Wg[i, j] <- 1; Wg[j, i] <- 1
  }
Lg <- diag(rowSums(Wg)) - Wg
blob <- exp(-((gidx$r - 4)^2 + (gidx$c - 4)^2) / 6); blob <- blob / sqrt(sum(blob^2))
chk  <- as.numeric(((gidx$r + gidx$c) %% 2) * 2 - 1); chk <- chk / sqrt(sum(chk^2))
Xg <- matrix(rnorm(nn), nn, 1) %*% t(chk) * 3 +
      matrix(rnorm(nn), nn, 1) %*% t(blob) * 15 +
      matrix(rnorm(nn * pp), nn, pp) * 0.5
f_sm  <- genpca(Xg, A = solve(diag(pp) + 2 * Lg), ncomp = 1, preproc = multivarious::center())
f_lap <- genpca(Xg, A = diag(pp) + 50 * Lg,      ncomp = 1, preproc = multivarious::center())
op <- par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
for (ff in list(list(f_sm, "Smoother metric: finds the blob"),
                list(f_lap, "Laplacian metric: finds the checkerboard"))) {
  v <- multivarious::components(ff[[1]])[, 1]
  image(matrix(v, gg, gg), axes = FALSE, main = ff[[2]], cex.main = 0.95,
        col = grey.colors(24, start = 0.95, end = 0.15))
}
par(op)

## ----ar1----------------------------------------------------------------------
rho     <- 0.7
n_t     <- 60
idx     <- 0:(n_t - 1)
Sigma_r <- outer(idx, idx, function(i, j) rho^abs(i - j))
M_ar1   <- solve(Sigma_r + 1e-3 * diag(n_t))

## ----rbf----------------------------------------------------------------------
coords <- as.matrix(expand.grid(x = 1:8, y = 1:8))
d2     <- as.matrix(dist(coords))^2
ell    <- 2
K      <- exp(-d2 / (2 * ell^2))
A_rbf  <- solve(K + 1e-3 * diag(nrow(K)))   # precision: emphasises fine scale
# A_smooth <- K + 1e-3 * diag(nrow(K))      # kernel itself: smooth loadings

## ----lap----------------------------------------------------------------------
W <- bandSparse(30, k = c(-1, 0, 1),
                diagonals = list(rep(0.2, 29),
                                 rep(1, 30),
                                 rep(0.2, 29)))
D     <- Diagonal(x = rowSums(W))
A_lap <- (D - W) + 1e-2 * Diagonal(nrow(W))
# A_smooth <- solve(A_lap)                  # smoother: spatially coherent loadings

## ----recipe-plots, echo = FALSE, fig.cap = "Three structured metrics, all shown in the precision (noise-whitening) orientation. Off-diagonal banding is what couples nearby rows or variables -- it tells GPCA 'treat these dimensions as related, not independent'.", fig.width = 8, fig.height = 3----
op <- par(mfrow = c(1, 3), mar = c(2, 2, 2, 1))
heat(M_ar1, "Inverse AR(1)")
heat(A_rbf, "Inverse RBF kernel")
heat(A_lap, "Regularised Laplacian")
par(op)

## ----mle-fit------------------------------------------------------------------
set.seed(1)
n_m <- 40; p_m <- 10
X_mle <- matrix(rnorm(n_m * p_m), n_m, p_m)
fit_mle <- gpca_mle(X_mle, ncomp = 2, max_iter = 6,
                    lambda = 1e-3, scale_fix = "none",
                    method = "eigen", verbose = FALSE)

metric_spectrum <- function(W) {
  ev <- eigen(as.matrix(W), symmetric = TRUE, only.values = TRUE)$values
  c(min = min(ev), max = max(ev), condition = max(ev) / min(ev))
}
signif(rbind(M = metric_spectrum(fit_mle$M),
             A = metric_spectrum(fit_mle$A)), 3)

## ----mle-plot, echo = FALSE, fig.cap = "Metric eigenvalues divided by their mean, on a log scale. The dashed line marks an identity-like spectrum; the learned row metric departs strongly from it.", fig.width = 7, fig.height = 3.5----
spectra <- lapply(list(M = fit_mle$M, A = fit_mle$A), function(W) {
  ev <- eigen(as.matrix(W), symmetric = TRUE, only.values = TRUE)$values
  ev / mean(ev)
})
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (nm in names(spectra)) {
  plot(spectra[[nm]], type = "b", pch = 19, log = "y",
       ylim = range(unlist(spectra)), xlab = "Eigenvalue index",
       ylab = "Eigenvalue / mean", main = paste("Learned", nm))
  abline(h = 1, lty = 2, col = "steelblue")
}
par(op)

## ----mle-spectrum-check, include = FALSE--------------------------------------
stopifnot(all(is.finite(unlist(spectra))), all(unlist(spectra) > 0))

## ----mle-progress-------------------------------------------------------------
data.frame(iteration = seq_along(fit_mle$loglik_path),
           penalized_loglik = fit_mle$loglik_path)

## ----spd-check, eval = FALSE--------------------------------------------------
# # Is the metric usable as-is, and what would a repair do to it?
# A_ok <- repair_metric(A, method = "ridge")
# attr(A_ok, "repair_report")
# 
# # Catch the repair warning programmatically inside a fit
# fit <- withCallingHandlers(
#   genpca(X, A = A, M = M, ncomp = 3, constraints_remedy = "ridge"),
#   genpca_metric_repaired = function(w) { print(w$report); invokeRestart("muffleWarning") }
# )

