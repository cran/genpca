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

grid_img <- function(z, g, main = "", zlim = range(z)) {
  image(matrix(z, g, g), axes = FALSE, main = main, cex.main = 0.95,
        zlim = zlim, col = colorRampPalette(c("steelblue", "white", "tomato"))(101))
  box(col = "grey80")
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

## ----simulate-----------------------------------------------------------------
set.seed(11)
g <- 16; p <- g * g; n <- 64
gr <- expand.grid(r = 1:g, c = 1:g)

# NOTE the orientation: spatial dimensions in ROWS, variables in COLUMNS
spat_cds <- rbind(gr$r, gr$c)
dim(spat_cds)

blob <- function(r0, c0, s = 1.6) {
  z <- exp(-((gr$r - r0)^2 + (gr$c - c0)^2) / (2 * s^2))
  z[z < 0.15] <- 0                       # compact support => sparse
  z / sqrt(sum(z^2))
}
v1 <- blob(5, 5); v2 <- blob(12, 12)

tt <- (0:(n - 1)) / n                    # orthogonal temporal profiles
u1 <- sin(2 * pi * tt); u1 <- u1 / sqrt(sum(u1^2))
u2 <- sin(4 * pi * tt); u2 <- u2 / sqrt(sum(u2^2))

signal <- 30 * tcrossprod(u1, v1) + 20 * tcrossprod(u2, v2)
X <- signal + matrix(rnorm(n * p, sd = 0.25), n, p)

c(true_support = sum(v1 != 0), of = p,
  SNR = round(norm(signal, "F") / norm(X - signal, "F"), 2))

## ----fit----------------------------------------------------------------------
fit <- sfpca(X, K = 2, spat_cds = spat_cds)
fit

## ----sparsity-----------------------------------------------------------------
V <- multivarious::components(fit)
c(nonzero_PC1 = sum(V[, 1] != 0), nonzero_PC2 = sum(V[, 2] != 0), of = p)

## ----compare------------------------------------------------------------------
pc <- prcomp(as.matrix(X), center = TRUE, rank. = 2)

rbind(
  sfpca = c(PC1 = abs(cor(V[, 1], v1)),        PC2 = abs(cor(V[, 2], v2))),
  pca   = c(PC1 = abs(cor(pc$rotation[, 1], v1)), PC2 = abs(cor(pc$rotation[, 2], v2)))
)

## ----recovery-plot, echo = FALSE, fig.cap = "Truth, sfpca, and PCA with signs aligned to truth and one common color scale: blue is negative, white is zero, red is positive. sfpca has a few extra nonzero sites in PC1; PCA has nonzero loadings throughout.", fig.width = 6, fig.height = 8.5----
# Align only the displayed vectors; retain the fitted factors for reconstruction.
align <- function(v, truth) if (sum(v * truth) < 0) -v else v
V_show <- cbind(align(V[, 1], v1), align(V[, 2], v2))
P_show <- cbind(align(pc$rotation[, 1], v1), align(pc$rotation[, 2], v2))
limit <- max(abs(c(v1, v2, V_show, P_show)))
op <- par(mfrow = c(3, 2), mar = c(1, 1, 2.5, 1))
grid_img(v1, g, "True pattern 1", c(-limit, limit))
grid_img(v2, g, "True pattern 2", c(-limit, limit))
grid_img(V_show[, 1], g, "sfpca loading 1", c(-limit, limit))
grid_img(V_show[, 2], g, "sfpca loading 2", c(-limit, limit))
grid_img(P_show[, 1], g, "PCA loading 1", c(-limit, limit))
grid_img(P_show[, 2], g, "PCA loading 2", c(-limit, limit))
par(op)

## ----temporal-----------------------------------------------------------------
U <- fit$ou
c(PC1 = abs(cor(U[, 1], u1)), PC2 = abs(cor(U[, 2], u2)))

## ----reconstruction-----------------------------------------------------------
Xhat  <- as.matrix(multivarious::reconstruct(fit))
pchat <- pc$x %*% t(pc$rotation) + matrix(pc$center, n, p, byrow = TRUE)

c(sfpca = norm(signal - Xhat,  "F") / norm(signal, "F"),
  pca   = norm(signal - pchat, "F") / norm(signal, "F"))

## ----recovery-check, include = FALSE------------------------------------------
stopifnot(all(is.finite(V_show)),
          cor(V_show[, 1], v1) > 0.95, cor(V_show[, 2], v2) > 0.95,
          norm(signal - Xhat, "F") < norm(signal - pchat, "F"))

## ----ablation-----------------------------------------------------------------
variants <- list(
  "defaults"                  = list(),
  "no sparsity (lambda_v = 0)" = list(lambda_v = 0),
  "no smoothing (alpha_v = 0)" = list(alpha_v = 0),
  "neither"                    = list(lambda_v = 0, alpha_v = 0),
  "heavy sparsity (lambda_v = 3)" = list(lambda_v = 3)
)

t(sapply(variants, function(extra) {
  f <- do.call(sfpca, c(list(X = X, K = 1, spat_cds = spat_cds), extra))
  v <- multivarious::components(f)[, 1]
  c(nonzero = sum(v != 0), cor_with_truth = round(abs(cor(v, v1)), 3))
}))

## ----selected-----------------------------------------------------------------
data.frame(
  component = 1:2,
  lambda_u  = signif(fit$lambda_u, 3), lambda_v = signif(fit$lambda_v, 3),
  alpha_u   = signif(fit$alpha_u, 3),  alpha_v  = signif(fit$alpha_v, 3)
)

## ----correlated-fit-----------------------------------------------------------
u2c <- sin(2 * pi * tt + 0.9); u2c <- u2c / sqrt(sum(u2c^2))
round(sum(u1 * u2c), 3)                 # the two profiles now overlap

set.seed(11)
Xc  <- 30 * tcrossprod(u1, v1) + 20 * tcrossprod(u2c, v2) +
       matrix(rnorm(n * p, sd = 0.25), n, p)
fc  <- sfpca(Xc, K = 2, spat_cds = spat_cds)
Uc  <- fc$ou; Vc <- multivarious::components(fc)

## ----orthogonality------------------------------------------------------------
round(crossprod(Uc), 3)     # would be the identity for a joint SVD
round(crossprod(Vc), 3)

## ----sdev-meaning-------------------------------------------------------------
Xm   <- as.matrix(Xc)
dc   <- multivarious::sdev(fc)
defl <- Xm - dc[1] * tcrossprod(Uc[, 1], Vc[, 1])

c(sdev_2      = dc[2],
  u2_X_v2     = as.numeric(t(Uc[, 2]) %*% Xm   %*% Vc[, 2]),   # does NOT match
  u2_Xdefl_v2 = as.numeric(t(Uc[, 2]) %*% defl %*% Vc[, 2]))   # matches

## ----orientation-trap, error = TRUE-------------------------------------------
try({
sfpca(X, K = 1, spat_cds = t(spat_cds))
})

