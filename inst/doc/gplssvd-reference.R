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

## ----quick-use----------------------------------------------------------------
set.seed(123)
N <- 150
shared <- rnorm(N)
X <- outer(shared, seq(0.5, 1.2, length.out = 8)) + matrix(rnorm(N * 8), N, 8)
Y <- outer(shared, seq(0.5, 1.2, length.out = 5)) + matrix(rnorm(N * 5), N, 5)
row_wt   <- diag(runif(N, 0.5, 1.5))
col_wt_x <- diag(runif(8, 0.8, 1.2))
col_wt_y <- diag(runif(5, 0.8, 1.2))

fit <- genpls(X, Y, ncomp = 2,
              preproc_x = multivarious::center(),
              preproc_y = multivarious::center(),
              Mx = row_wt, My = row_wt,
              Ax = col_wt_x, Ay = col_wt_y)
round(fit$d, 3)

## ----quick-project------------------------------------------------------------
Sx <- multivarious::project(fit, X)
Sy <- multivarious::project(fit, Y, source = "Y")
cor(Sx[, 1], Sy[, 1])

## ----quick-plot, echo = FALSE, fig.cap = "Leading projected coordinates of the two blocks, which share a simulated signal.", fig.height = 3.5----
plot(Sx[, 1], Sy[, 1], pch = 19, col = "steelblue",
     xlab = "X latent coordinate 1", ylab = "Y latent coordinate 1")

## ----quick-project-check, include = FALSE-------------------------------------
stopifnot(identical(dim(Sx), c(150L, 2L)), identical(dim(Sy), c(150L, 2L)),
          all(is.finite(Sx)), all(is.finite(Sy)),
          abs(cor(Sx[, 1], Sy[, 1])) > 0.6)

## ----ref-sqrt-----------------------------------------------------------------
psd_sqrt <- function(W, n) {
  if (is.null(W)) return(list(h = diag(n), hinv = diag(n), full = diag(n)))
  W    <- as.matrix(W)
  e    <- eigen(W, symmetric = TRUE)
  lam  <- pmax(e$values, 0)                    # PSD assumed; clip numerical negatives
  half <- function(f) e$vectors %*% (f * t(e$vectors))
  list(h    = half(sqrt(lam)),                 # W^{1/2}
       hinv = half(ifelse(lam > 0, 1 / sqrt(lam), 0)),  # W^{-1/2}, pseudo
       full = W)
}

## ----ref-impl-----------------------------------------------------------------
dense_gplssvd_ref <- function(X, Y, MX = NULL, MY = NULL,
                              WX = NULL, WY = NULL, k = NULL,
                              center = FALSE, scale = FALSE) {
  X <- scale(as.matrix(X), center = center, scale = scale)
  Y <- scale(as.matrix(Y), center = center, scale = scale)
  stopifnot(nrow(X) == nrow(Y))

  mx <- psd_sqrt(MX, nrow(X)); wx <- psd_sqrt(WX, ncol(X))
  my <- psd_sqrt(MY, nrow(Y)); wy <- psd_sqrt(WY, ncol(Y))

  # whiten both blocks, then SVD their cross-product
  S  <- crossprod(mx$h %*% X %*% wx$h, my$h %*% Y %*% wy$h)
  sv <- svd(S)
  keep <- seq_len(if (is.null(k)) length(sv$d) else min(k, length(sv$d)))

  # unwhiten: saliences are W^{-1/2} u, so that p' WX p = I
  p <- wx$hinv %*% sv$u[, keep, drop = FALSE]
  q <- wy$hinv %*% sv$v[, keep, drop = FALSE]
  D <- diag(sv$d[keep], nrow = length(keep))

  list(d  = sv$d[keep], p = p, q = q,
       fi = wx$full %*% p %*% D,             # factor scores
       fj = wy$full %*% q %*% D,
       lx = mx$h %*% X %*% wx$full %*% p,    # latent variables
       ly = my$h %*% Y %*% wy$full %*% q)
}

## ----example------------------------------------------------------------------
set.seed(1)
N <- 20; I <- 8; J <- 6
X  <- matrix(rnorm(N * I), N, I)
Y  <- matrix(rnorm(N * J), N, J)
MX <- diag(runif(N, .5, 1.5))
MY <- diag(runif(N, .5, 1.5))
WX <- diag(runif(I, .5, 1.5))
WY <- diag(runif(J, .5, 1.5))

ref <- dense_gplssvd_ref(X, Y, MX, MY, WX, WY,
                         k = 3, center = TRUE, scale = FALSE)
op  <- gplssvd_op(X, Y,
                  XLW = MX, YLW = MY,
                  XRW = WX, YRW = WY,
                  k = 3, center = TRUE, scale = FALSE)

all.equal(ref$d, op$d, tolerance = 1e-6)
all.equal(diag(crossprod(op$lx, op$ly)), op$d, tolerance = 1e-6)
round(op$d, 4)

## ----ref-vs-op-plot, echo = FALSE, fig.cap = "Reference vs operator singular values agree to plotting precision (left). Latent variables show the expected diagonal cross-product structure (right).", fig.width = 7, fig.height = 3.5----
op_par <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
plot(ref$d, type = "b", pch = 19, col = "grey30",
     xlab = "Component", ylab = "Singular value",
     main = "Singular values")
lines(op$d, type = "b", pch = 21, col = "steelblue")
legend("topright", legend = c("reference", "operator"),
       col = c("grey30", "steelblue"), pch = c(19, 21),
       bty = "n", cex = 0.85)

cp <- as.matrix(crossprod(op$lx, op$ly))
image(t(cp)[, ncol(cp):1], axes = FALSE,
      main = "t(Lx) %*% Ly",
      col = grey.colors(20, start = 0.95, end = 0.2))
par(op_par)

