context("dual-path (n < p) correctness with non-identity metrics")

library(Matrix)
library(multivarious)

# Regression tests for a bug where the eigen method's dual path (n < p)
# omitted the row metric when recovering the second factor, producing
# corrupted ou/scores/reconstruction whenever A != I even though singular
# values were correct. These tests assert the invariants that the old
# differential tests missed: metric-orthonormality of ou/ov for EVERY
# backend, and full-rank reconstruction of the input.

make_sparse_spd <- function(m, density = 0.05, seed = 1) {
  set.seed(seed)
  B <- rsparsematrix(m, m, density = density)
  forceSymmetric(crossprod(B) + Diagonal(m))
}

expect_gmd_invariants <- function(fit, X, M, A, tol_orth = 1e-6, tol_recon = 1e-6,
                                  full_rank = FALSE) {
  ou <- as.matrix(fit$ou)
  ov <- as.matrix(fit$ov)
  k <- length(fit$sdev)
  expect_gt(k, 0)
  expect_lt(max(abs(as.matrix(crossprod(ou, M %*% ou)) - diag(k))), tol_orth)
  expect_lt(max(abs(as.matrix(crossprod(ov, A %*% ov)) - diag(k))), tol_orth)
  if (full_rank) {
    Xrec <- ou %*% (fit$sdev * t(ov))
    expect_lt(max(abs(Xrec - as.matrix(X))), tol_recon)
  }
}

test_that("eigen dual path (n < p, A != I diagonal) satisfies GMD invariants", {
  set.seed(10)
  n <- 8; p <- 20
  X <- matrix(rnorm(n * p), n, p)
  A <- Diagonal(x = runif(p, 0.5, 3))
  M <- Diagonal(x = runif(n, 0.5, 2))

  fit <- genpca(X, A = A, M = M, ncomp = n, method = "eigen",
                preproc = multivarious::pass())
  expect_gmd_invariants(fit, X, M, A, full_rank = TRUE, tol_recon = 1e-8)
})

test_that("eigen dual path (n < p, general sparse SPD metrics) satisfies GMD invariants", {
  set.seed(11)
  n <- 12; p <- 40
  X <- matrix(rnorm(n * p), n, p)
  M <- make_sparse_spd(n, density = 0.2, seed = 2)
  A <- make_sparse_spd(p, density = 0.1, seed = 3)

  fit <- genpca(X, A = A, M = M, ncomp = n, method = "eigen",
                preproc = multivarious::pass())
  expect_gmd_invariants(fit, X, M, A, full_rank = TRUE, tol_recon = 1e-8)
})

test_that("all backends agree on the wide (n < p, A != I) problem", {
  set.seed(12)
  n <- 15; p <- 60; k <- 5
  X <- matrix(rnorm(n * p), n, p)
  A <- Diagonal(x = runif(p, 0.5, 3))
  M <- Diagonal(x = runif(n, 0.5, 2))

  fit_eig <- genpca(X, A = A, M = M, ncomp = k, method = "eigen",
                    preproc = multivarious::pass())
  fit_spc <- genpca(X, A = A, M = M, ncomp = k, method = "spectra",
                    preproc = multivarious::pass())

  expect_equal(fit_eig$sdev, fit_spc$sdev, tolerance = 1e-6)

  # scores must agree up to sign, not just sdev
  sc_e <- as.matrix(scores(fit_eig))
  sc_s <- as.matrix(scores(fit_spc))
  sgn <- sign(colSums(sc_e * sc_s))
  expect_lt(max(abs(sweep(sc_s, 2, sgn, `*`) - sc_e)), 1e-5)

  # and both satisfy the metric-orthonormality invariants
  expect_gmd_invariants(fit_eig, X, M, A)
  expect_gmd_invariants(fit_spc, X, M, A)
})

test_that("eigen method returns all min(n, p) components at full rank", {
  set.seed(13)
  # tall case: full rank = p
  Xt <- matrix(rnorm(30 * 10), 30, 10)
  fit_t <- genpca(Xt, ncomp = 10, method = "eigen", preproc = multivarious::pass())
  expect_equal(length(fit_t$sdev), 10)

  # wide case: full rank = n
  Xw <- matrix(rnorm(10 * 30), 10, 30)
  fit_w <- genpca(Xw, ncomp = 10, method = "eigen", preproc = multivarious::pass())
  expect_equal(length(fit_w$sdev), 10)
})

test_that("C++ deflation handles a zero matrix consistently with the R path", {
  fit <- suppressWarnings(
    genpca(matrix(0, 6, 4), ncomp = 3, method = "deflation", use_cpp = TRUE)
  )
  expect_equal(multivarious::ncomp(fit), 0)
  expect_equal(dim(multivarious::scores(fit)), c(6, 0))
})

test_that("genpca_cov_gmd diagonal fast path matches dense-eigen reference", {
  set.seed(14)
  p <- 25
  Z <- matrix(rnorm(100 * p), 100, p)
  C <- crossprod(Z) / 100
  w <- runif(p, 0.5, 2)

  res <- genpca:::genpca_cov_gmd(C, R = w, ncomp = 5)

  # reference: eigenvalues of D^{1/2} C D^{1/2}
  B <- outer(sqrt(w), sqrt(w)) * as.matrix(C)
  lam_ref <- sort(eigen(B, symmetric = TRUE, only.values = TRUE)$values,
                  decreasing = TRUE)[1:5]
  expect_equal(as.numeric(res$lambda), lam_ref, tolerance = 1e-8)

  # loadings R-orthonormal
  V <- as.matrix(res$v)
  expect_lt(max(abs(crossprod(V, w * V) - diag(5))), 1e-8)
})
