library(testthat)
library(genpca)

make_chain_precision <- function(m, rho = 0.32, jitter = 0.7) {
  T <- matrix(0, m, m)
  diag(T) <- 1 + jitter
  for (i in seq_len(m - 1L)) {
    T[i, i + 1L] <- -rho
    T[i + 1L, i] <- -rho
  }
  T
}

make_block_precision <- function(p, b = 4, in_block = 0.24, jitter = 0.6) {
  T <- matrix(0, p, p)
  nb <- ceiling(p / b)
  for (k in seq_len(nb)) {
    idx <- ((k - 1L) * b + 1L):min(k * b, p)
    T[idx, idx] <- in_block
  }
  diag(T) <- diag(T) + 1 + jitter
  T
}

matnorm_noise <- function(n, p, Omega, Sigma) {
  Z <- matrix(rnorm(n * p), n, p)
  chol(Omega) %*% Z %*% t(chol(Sigma))
}

offdiag_support_metrics <- function(Theta_true, Theta_est, thr = 1e-3) {
  mask <- row(Theta_true) != col(Theta_true)
  truth <- abs(Theta_true[mask]) > 1e-9
  pred <- abs(Theta_est[mask]) > thr
  tp <- sum(truth & pred)
  fn <- sum(truth & !pred)
  fp <- sum(!truth & pred)
  tn <- sum(!truth & !pred)
  c(
    tpr = tp / max(1, tp + fn),
    fpr = fp / max(1, fp + tn)
  )
}

offdiag_mean <- function(M) {
  A <- as.matrix(M)
  mean(abs(A[row(A) != col(A)]))
}

test_that("glasso row update recovers chain precision from known row structure", {
  set.seed(1003)
  n <- 22
  p <- 280

  Theta_row_true <- make_chain_precision(n)
  E <- matnorm_noise(n, p, Omega = solve(Theta_row_true), Sigma = diag(p))
  S1 <- (E %*% t(E)) / p

  row_fit <- genpca:::.mnpca_glasso_admm(
    S = S1,
    lambda = 0.08,
    maxit = 400,
    tol = 1e-4,
    block_screen = TRUE,
    penalize_diagonal = FALSE
  )

  m <- offdiag_support_metrics(Theta_row_true, row_fit$Theta, thr = 1e-3)
  expect_gte(m[["tpr"]], 0.75)
  expect_lte(m[["fpr"]], 0.10)
  expect_true(genpca:::is_spd(as.matrix(row_fit$Theta)))
})

test_that("glasso column update recovers block precision from known column structure", {
  set.seed(2005)
  n <- 280
  p <- 22

  Theta_col_true <- make_block_precision(p, b = 4)
  E <- matnorm_noise(n, p, Omega = diag(n), Sigma = solve(Theta_col_true))
  S2 <- (t(E) %*% E) / n

  col_fit <- genpca:::.mnpca_glasso_admm(
    S = S2,
    lambda = 0.05,
    maxit = 400,
    tol = 1e-4,
    block_screen = TRUE,
    penalize_diagonal = FALSE
  )

  m <- offdiag_support_metrics(Theta_col_true, col_fit$Theta, thr = 1e-3)
  expect_gte(m[["tpr"]], 0.50)
  expect_lte(m[["fpr"]], 0.20)
  expect_true(genpca:::is_spd(as.matrix(col_fit$Theta)))
})

test_that("end-to-end mnpca_mrl distinguishes row-only structured noise", {
  set.seed(3004)
  n <- 30
  p <- 24
  r <- 3

  Theta_row_true <- make_chain_precision(n)
  signal <- matrix(rnorm(n * r), n, r) %*% t(matrix(rnorm(p * r), p, r))
  Y <- signal + matnorm_noise(n, p, Omega = solve(Theta_row_true), Sigma = diag(p))

  fit <- mnpca_mrl(
    Y,
    ncomp = r,
    center = FALSE,
    update_precisions = TRUE,
    lambda_row = 0.08,
    lambda_col = 0.12,
    max_outer = 8,
    max_inner = 7,
    gl_maxit = 140,
    gl_tol = 1e-3,
    tol = 1e-4,
    as_sparse_precision = FALSE,
    verbose = FALSE
  )

  row_off <- offdiag_mean(fit$Theta_row)
  col_off <- offdiag_mean(fit$Theta_col)

  expect_gt(row_off, 0.01)
  expect_lt(col_off, 0.01)
  expect_gt(row_off, 2.5 * col_off)
})

test_that("end-to-end mnpca_mrl distinguishes column-only structured noise", {
  set.seed(3506)
  n <- 30
  p <- 24
  r <- 3

  Theta_col_true <- make_block_precision(p, b = 4)
  signal <- matrix(rnorm(n * r), n, r) %*% t(matrix(rnorm(p * r), p, r))
  Y <- signal + matnorm_noise(n, p, Omega = diag(n), Sigma = solve(Theta_col_true))

  fit <- mnpca_mrl(
    Y,
    ncomp = r,
    center = FALSE,
    update_precisions = TRUE,
    lambda_row = 0.12,
    lambda_col = 0.08,
    max_outer = 8,
    max_inner = 7,
    gl_maxit = 140,
    gl_tol = 1e-3,
    tol = 1e-4,
    as_sparse_precision = FALSE,
    verbose = FALSE
  )

  row_off <- offdiag_mean(fit$Theta_row)
  col_off <- offdiag_mean(fit$Theta_col)

  expect_lt(row_off, 0.01)
  expect_gt(col_off, 0.007)
  expect_gt(col_off, 2 * row_off)
})
