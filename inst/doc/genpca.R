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

## ----quick-start--------------------------------------------------------------
data("USArrests")
X <- as.matrix(USArrests[, c("Murder", "Assault", "Rape")])
col_sd <- apply(X, 2, sd)
A <- Diagonal(x = 1 / col_sd^2)
fit <- genpca(X, A = A, ncomp = 2,
              preproc = multivarious::center())

S <- multivarious::scores(fit)       # 50 states x 2 components
V <- multivarious::components(fit)   # 3 variables x 2 projection weights
round(cor(X, S), 2)

## ----quick-plot, echo = FALSE, fig.cap = "States in the two-component space after inverse-variance weighting of the crime variables. Six states with extreme scores are labelled.", fig.height = 4.5----
plot(S, pch = 19, col = "grey50",
     xlab = "PC1", ylab = "PC2", main = "Standardized crime variables",
     xlim = extendrange(S[, 1], f = 0.15),
     ylim = extendrange(S[, 2], f = 0.15))
label <- unique(c(order(S[, 1])[c(1, 2, 49, 50)], which.min(S[, 2]), which.max(S[, 2])))
text(S[label, 1], S[label, 2], labels = state.abb[label], pos = 3, cex = 0.8)

## ----scaling-comparison-------------------------------------------------------
fit_raw <- genpca(X, ncomp = 2, preproc = multivarious::center())
round(cbind(
  unscaled = abs(cor(X, multivarious::scores(fit_raw)[, 1])),
  scaled   = abs(cor(X, S[, 1]))
), 2)

## ----quick-check, include = FALSE---------------------------------------------
stopifnot(identical(dim(S), c(50L, 2L)), all(is.finite(S)),
          all(abs(cor(X, S[, 1])) > 0.7))

## ----object-flow--------------------------------------------------------------
Xhat <- multivarious::reconstruct(fit)   # 50 x 3, back in original units
scores_again <- multivarious::project(fit, X)
max(abs(S - scores_again))              # numerical zero on training rows

## ----row-weights--------------------------------------------------------------
urban_wt <- USArrests$UrbanPop / mean(USArrests$UrbanPop)
M <- Diagonal(x = urban_wt)
fit_row <- genpca(X, M = M, A = A, ncomp = 2,
                  preproc = multivarious::center())
round(cor(X, multivarious::scores(fit_row)), 2)

## ----scree, fig.cap = "Fraction of total weighted variation for all three available components.", fig.height = 3----
fit_all <- genpca(X, A = A, ncomp = ncol(X),
                  preproc = multivarious::center())
share <- fit_all$sdev^2 / sum(fit_all$sdev^2)
barplot(share, names.arg = paste0("PC", seq_along(share)),
        ylim = c(0, 1), ylab = "Fraction of weighted variation",
        col = "grey60", border = NA)

