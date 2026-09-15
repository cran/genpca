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

## ----graph-setup--------------------------------------------------------------
g <- 12; p <- g * g; n <- 150
gi <- expand.grid(r = 1:g, c = 1:g)
W <- matrix(0, p, p)
for (i in 1:p) for (j in 1:p)
  if (i < j && abs(gi$r[i] - gi$r[j]) + abs(gi$c[i] - gi$c[j]) == 1) {
    W[i, j] <- 1; W[j, i] <- 1
  }
L  <- diag(rowSums(W)) - W
eL <- eigen(L, symmetric = TRUE)
Q  <- eL$vectors
lam <- pmax(eL$values, 0)

# build a metric from a transfer function of the Laplacian
spec <- function(f) Q %*% (f(lam) * t(Q))

## ----profiles, echo = FALSE, fig.cap = "Four transfer functions on the same graph, each normalized to a maximum of one. Lower Laplacian eigenvalues represent smoother patterns; there is no single boundary between smooth and rough.", fig.height = 3.6----
ls <- seq(0, max(lam), length.out = 200)
prof <- cbind(smoother  = 1 / (1 + 6 * ls),
              precision = 1 / (0.02 + 1 / (1 + 6 * ls)),
              unbounded = 1 + 6 * ls,
              bandpass  = exp(-((ls - 3.5)^2) / 1.5))
prof <- sweep(prof, 2, apply(prof, 2, max), "/")
matplot(ls, prof, type = "l", lty = 1, lwd = 2,
        col = c("steelblue", "tomato", "grey30", "darkolivegreen"),
        xlab = expression(paste("Laplacian eigenvalue  ", lambda,
                                "   (smooth ", symbol("\256"), " rough)")),
        ylab = "relative weight f(lambda)")
legend("right", c("smoother", "bounded precision", "unbounded I + aL", "band-pass"),
       col = c("steelblue", "tomato", "grey30", "darkolivegreen"),
       lty = 1, lwd = 2, bty = "n", cex = 0.8)

## ----adjoin, eval = requireNamespace("adjoin", quietly = TRUE)----------------
cds  <- as.matrix(expand.grid(x = 1:8, y = 1:8))
Aadj <- adjoin::spatial_adjacency(cds, nnk = 8, weight_mode = "heat", sigma = 1.5)
Alap <- adjoin::spatial_laplacian(cds, nnk = 8, weight_mode = "heat", sigma = 1.5)

# check the two traps before using either as a metric
range(eigen(as.matrix(Aadj), symmetric = TRUE, only.values = TRUE)$values)
range(eigen(as.matrix(Alap), symmetric = TRUE, only.values = TRUE)$values)

## ----twobytwo-----------------------------------------------------------------
nrm  <- function(v) v / sqrt(sum(v^2))
# how much of the planted pattern is captured by the fitted subspace?
recov <- function(V, truth) {
  V <- qr.Q(qr(as.matrix(V)))
  sqrt(sum((t(V) %*% truth)^2)) / sqrt(sum(truth^2))
}

smooth_pat <- nrm(exp(-((gi$r - 4)^2 + (gi$c - 4)^2) / 6))     # a blob
fine_pat   <- nrm(Q[, which.min(abs(lam - median(lam)))])      # a mid-frequency mode

smooth_noise <- function(n) {
  Z <- matrix(rnorm(n * p), n, p) %*% spec(function(l) sqrt(1 / (1 + 4 * l)))
  Z / sqrt(mean(Z^2))
}
fine_noise <- function(n) {
  Z <- matrix(rnorm(n * p), n, p) %*% spec(function(l) sqrt((l + .5) / max(lam)))
  Z / sqrt(mean(Z^2))
}

A_smoother  <- spec(function(l) 1 / (1 + 6 * l))
A_precision <- spec(function(l) 1 / (0.02 + 1 / (1 + 6 * l)))

set.seed(909)
reps <- 8
grid <- expand.grid(signal = c("smooth", "fine"), noise = c("smooth", "fine"),
                    stringsAsFactors = FALSE)
out <- t(apply(grid, 1, function(row) {
  pat <- if (row[["signal"]] == "smooth") smooth_pat else fine_pat
  acc <- c(0, 0, 0)
  for (r in seq_len(reps)) {
    E <- if (row[["noise"]] == "smooth") smooth_noise(n) else fine_noise(n)
    X <- scale(matrix(rnorm(n), n, 1) %*% t(pat) * 1.1 + E, scale = FALSE)
    for (k in 1:3) {
      A <- list(diag(p), A_smoother, A_precision)[[k]]
      fit <- genpca(X, A = A, ncomp = 1, preproc = multivarious::pass())
      acc[k] <- acc[k] + recov(fit$ov, pat) / reps
    }
  }
  acc
}))
dimnames(out) <- list(paste(grid$signal, "signal /", grid$noise, "noise"),
                      c("identity", "smoother", "precision"))
round(out, 3)

## ----bounded------------------------------------------------------------------
A_unbounded <- spec(function(l) 1 + 6 * l)

set.seed(78)
acc <- c(0, 0, 0); reps <- 8
for (r in seq_len(reps)) {
  X <- scale(matrix(rnorm(n), n, 1) %*% t(fine_pat) * 1.1 +
             smooth_noise(n) * 0.8 +
             matrix(rnorm(n * p), n, p) * 0.8, scale = FALSE)   # + thermal
  for (k in 1:3) {
    A <- list(diag(p), A_precision, A_unbounded)[[k]]
    fit <- genpca(X, A = A, ncomp = 1, preproc = multivarious::pass())
    acc[k] <- acc[k] + recov(fit$ov, fine_pat) / reps
  }
}
setNames(round(acc, 3), c("identity", "bounded precision", "unbounded I + 6L"))

## ----rowmetric----------------------------------------------------------------
rho   <- 0.85
Sig_t <- outer(0:(n - 1), 0:(n - 1), function(i, j) rho^abs(i - j))
M_ar  <- solve(Sig_t + 1e-6 * diag(n))
task  <- scale(sin(2 * pi * (1:n) / 7))

set.seed(303)
acc <- c(0, 0); reps <- 8
for (r in seq_len(reps)) {
  E <- t(chol(Sig_t)) %*% matrix(rnorm(n * p), n, p)      # AR(1) in time
  X <- scale(task %*% t(smooth_pat) + E, scale = FALSE)
  acc[1] <- acc[1] + recov(genpca(X, ncomp = 1,
                                  preproc = multivarious::pass())$ov, smooth_pat) / reps
  acc[2] <- acc[2] + recov(genpca(X, M = M_ar, ncomp = 1,
                                  preproc = multivarious::pass())$ov, smooth_pat) / reps
}
setNames(round(acc, 3), c("no row metric", "AR(1) precision M"))

