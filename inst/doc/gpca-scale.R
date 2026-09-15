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

## ----backend-bench------------------------------------------------------------
set.seed(11)
n <- 150; p <- 60
X <- matrix(rnorm(n * p), n, p)

t_eig <- system.time(
  fit_eig <- genpca(X, ncomp = 8, method = "eigen",
                    preproc = multivarious::center())
)
t_rnd <- system.time(
  fit_rnd <- genpca(X, ncomp = 8, method = "randomized",
                    preproc = multivarious::center())
)
data.frame(method = c("eigen", "randomized"),
           elapsed = c(t_eig["elapsed"], t_rnd["elapsed"]),
           top_sv  = c(fit_eig$sdev[1], fit_rnd$sdev[1]),
           max_relative_error = c(0, max(abs(fit_rnd$sdev / fit_eig$sdev - 1))))

## ----backend-plot, echo = FALSE, fig.cap = "The randomized approximation underestimates the reference singular values on this full-rank example. The table reports the largest relative difference.", fig.height = 3.5----
plot(fit_eig$sdev, type = "b", pch = 19, col = "grey30",
     ylim = range(c(fit_eig$sdev, fit_rnd$sdev)),
     xlab = "Component", ylab = "Singular value",
     main = "Backend comparison")
lines(fit_rnd$sdev, type = "b", pch = 21, col = "steelblue")
legend("topright", legend = c("eigen", "randomized"),
       col = c("grey30", "steelblue"), pch = c(19, 21),
       bty = "n", cex = 0.85)

## ----spectra------------------------------------------------------------------
set.seed(42)
n <- 300; p <- 200
X_sparse <- rsparsematrix(n, p, density = 0.01)

# Sparse tridiagonal row/column metrics (mild AR(1)-style coupling)
M_sp <- bandSparse(n, k = c(-1, 0, 1),
                   diagonals = list(rep(0.1, n - 1), rep(1, n), rep(0.1, n - 1)))
A_sp <- bandSparse(p, k = c(-1, 0, 1),
                   diagonals = list(rep(0.1, p - 1), rep(1, p), rep(0.1, p - 1)))

fit_sp <- genpca(X_sparse, M = M_sp, A = A_sp, ncomp = 5, method = "spectra",
                 preproc = multivarious::pass())
fit_sp$sdev

## ----cov----------------------------------------------------------------------
set.seed(123)
n <- 100; p <- 15
X <- matrix(rnorm(n * p), n, p)
M <- diag(runif(n, 0.8, 1.2))
A <- diag(runif(p, 0.7, 1.3))
C <- t(X) %*% M %*% X
fit_cov <- genpca_cov(C, R = A, ncomp = 5, method = "gmd")
fit_cov$d

## ----cov-plot, echo = FALSE, fig.cap = "Singular values from the covariance-only fit.", fig.height = 3----
barplot(fit_cov$d, names.arg = paste0("PC", seq_along(fit_cov$d)),
        col = "grey60", border = NA, ylab = "Singular value")

## ----oos----------------------------------------------------------------------
set.seed(7)
X <- matrix(rnorm(200 * 30), 200, 30)
fit <- genpca(X[1:150, ], ncomp = 4,
              preproc = multivarious::center())
scores_test <- multivarious::project(fit, X[151:200, ])
head(scores_test, 4)

## ----oos-plot, echo = FALSE, fig.cap = "Training scores (grey) and out-of-sample scores (blue) projected into the same component space.", fig.height = 4----
S_train <- multivarious::scores(fit)
plot(rbind(S_train, scores_test)[, 1:2], type = "n",
     xlab = "PC1", ylab = "PC2",
     main = "Training vs out-of-sample")
points(S_train[, 1], S_train[, 2], pch = 19, col = "grey60")
points(scores_test[, 1], scores_test[, 2], pch = 19, col = "steelblue")
legend("topright", legend = c("Train", "OOS"),
       col = c("grey60", "steelblue"), pch = 19,
       bty = "n", cex = 0.85)

