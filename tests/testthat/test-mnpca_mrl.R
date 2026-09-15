library(testthat)
library(genpca)

test_that("mnpca_mrl returns finite outputs with SPD precisions", {
  set.seed(11)
  Y <- matrix(rnorm(80), 10, 8)

  fit <- mnpca_mrl(
    Y,
    ncomp = 3,
    lambda_row = 0.05,
    lambda_col = 0.05,
    max_outer = 4,
    max_inner = 4,
    gl_maxit = 80,
    gl_tol = 1e-3,
    tol = 1e-4,
    verbose = FALSE
  )

  expect_equal(dim(fit$X), c(nrow(Y), 3))
  expect_equal(dim(fit$W), c(ncol(Y), 3))
  expect_equal(dim(fit$Theta_row), c(nrow(Y), nrow(Y)))
  expect_equal(dim(fit$Theta_col), c(ncol(Y), ncol(Y)))
  expect_equal(dim(fit$fitted), dim(Y))
  expect_true(all(is.finite(fit$fitted_centered)))
  expect_true(all(is.finite(fit$residual_centered)))
  expect_true(length(fit$objective_path) >= 1)
  expect_true(all(is.finite(fit$objective_path)))
  expect_true(genpca:::is_spd(as.matrix(fit$Theta_row)))
  expect_true(genpca:::is_spd(as.matrix(fit$Theta_col)))
})

test_that("residual-free scatter updates match explicit residual formulas", {
  set.seed(22)
  n <- 7
  p <- 6
  r <- 2

  Y <- matrix(rnorm(n * p), n, p)
  X <- matrix(rnorm(n * r), n, r)
  W <- matrix(rnorm(p * r), p, r)

  Ar <- crossprod(matrix(rnorm(n * n), n, n)) / n + diag(0.5, n)
  Ac <- crossprod(matrix(rnorm(p * p), p, p)) / p + diag(0.5, p)

  fast <- genpca:::.mnpca_scatter_no_residual(
    Y = Y,
    X = X,
    W = W,
    Theta_row = Ar,
    Theta_col = Ac,
    jitter = 0
  )

  E <- Y - X %*% t(W)
  S1_ref <- (E %*% Ac %*% t(E)) / p
  S2_ref <- (t(E) %*% Ar %*% E) / n
  S1_ref <- 0.5 * (S1_ref + t(S1_ref))
  S2_ref <- 0.5 * (S2_ref + t(S2_ref))

  expect_equal(fast$S1, S1_ref, tolerance = 1e-8)
  expect_equal(fast$S2, S2_ref, tolerance = 1e-8)
})

test_that("mnpca_mrl ALS-only mode matches truncated SVD reconstruction", {
  set.seed(33)
  Y <- matrix(rnorm(12 * 9), 12, 9)
  r <- 3

  fit <- mnpca_mrl(
    Y,
    ncomp = r,
    center = FALSE,
    update_precisions = FALSE,
    lambda_row = 0,
    lambda_col = 0,
    max_outer = 3,
    max_inner = 60,
    tol = 1e-8,
    as_sparse_precision = FALSE,
    verbose = FALSE
  )

  sv <- svd(Y)
  Y_opt <- sv$u[, seq_len(r), drop = FALSE] %*%
    diag(sv$d[seq_len(r)], nrow = r, ncol = r) %*%
    t(sv$v[, seq_len(r), drop = FALSE])

  err_fit <- norm(Y - fit$fitted_centered, type = "F")
  err_opt <- norm(Y - Y_opt, type = "F")

  expect_lte(err_fit, err_opt * (1 + 1e-3) + 1e-8)
})

test_that("larger penalties shrink precision off-diagonals", {
  set.seed(44)
  Y <- matrix(rnorm(90), 10, 9)

  fit_low <- mnpca_mrl(
    Y,
    ncomp = 2,
    lambda_row = 0.01,
    lambda_col = 0.01,
    max_outer = 3,
    max_inner = 4,
    gl_maxit = 80,
    gl_tol = 1e-3,
    tol = 1e-4,
    verbose = FALSE
  )

  fit_high <- mnpca_mrl(
    Y,
    ncomp = 2,
    lambda_row = 0.25,
    lambda_col = 0.25,
    max_outer = 3,
    max_inner = 4,
    gl_maxit = 80,
    gl_tol = 1e-3,
    tol = 1e-4,
    verbose = FALSE
  )

  offdiag_mean <- function(M) {
    A <- as.matrix(M)
    A[row(A) == col(A)] <- NA_real_
    mean(abs(A), na.rm = TRUE)
  }

  low_mean <- offdiag_mean(fit_low$Theta_row) + offdiag_mean(fit_low$Theta_col)
  high_mean <- offdiag_mean(fit_high$Theta_row) + offdiag_mean(fit_high$Theta_col)

  expect_lte(high_mean, low_mean + 1e-8)
})

test_that("mnpca_mrl matches genpca in identity-metric ALS-only regime", {
  set.seed(55)
  n <- 40
  p <- 25
  r <- 5
  Y <- matrix(rnorm(n * p), n, p)

  fit_gpca <- genpca(
    Y,
    ncomp = r,
    preproc = multivarious::center(),
    method = "eigen"
  )
  Yhat_gpca <- suppressWarnings(multivarious::reconstruct(fit_gpca))

  fit_mn <- mnpca_mrl(
    Y,
    ncomp = r,
    center = TRUE,
    update_precisions = FALSE,
    lambda_row = 0,
    lambda_col = 0,
    max_outer = 5,
    max_inner = 120,
    tol = 1e-8,
    as_sparse_precision = FALSE,
    verbose = FALSE
  )
  Yhat_mn <- fit_mn$fitted

  rel_recon <- norm(Yhat_gpca - Yhat_mn, type = "F") / max(1, norm(Yhat_gpca, type = "F"))
  expect_lte(rel_recon, 1e-6)

  Yc <- scale(Y, center = TRUE, scale = FALSE)
  err_gpca <- norm(Yc - scale(Yhat_gpca, center = TRUE, scale = FALSE), type = "F")
  err_mn <- norm(Yc - scale(Yhat_mn, center = TRUE, scale = FALSE), type = "F")
  expect_equal(err_mn, err_gpca, tolerance = 1e-6)

  Vg <- as.matrix(fit_gpca$ov[, seq_len(r), drop = FALSE])
  Vm <- qr.Q(qr(fit_mn$W))
  cosines <- svd(t(Vg) %*% Vm)$d
  expect_gte(min(cosines), 0.9999)
})
