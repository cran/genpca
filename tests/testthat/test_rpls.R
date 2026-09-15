# File: tests/testthat/test-rpls.R

library(testthat)
library(multivarious)

test_that("rpls validates input dimensions and lambda", {
  set.seed(1)
  X <- matrix(rnorm(10 * 3), 10, 3)
  Y_bad <- matrix(rnorm(11 * 2), 11, 2)
  expect_error(rpls(X, Y_bad), "same number of rows")

  Y_good <- matrix(rnorm(10 * 2), 10, 2)
  expect_error(rpls(X, Y_good, lambda = -0.1), "non-negative")
  expect_error(rpls(X, Y_good, lambda = "a"), "numeric")
})

test_that("rpls with default arguments works correctly for a small example", {
  # Generate synthetic data
  set.seed(42)
  n <- 20
  p <- 6
  q <- 3
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * q), nrow = n, ncol = q)

  # Center columns (typical for PLS)
  X <- scale(X, center = TRUE, scale = FALSE)
  Y <- scale(Y, center = TRUE, scale = FALSE)

  # Fit RPLS with default arguments
  fit <- rpls(X, Y, K = 2, lambda = 0.1)

  # Check that returned object has the correct classes
  expect_s3_class(fit, c("rpls", "cross_projector", "projector"))

  # Check that the number of components matches K (or less if it broke early)
  expect_true(ncol(fit$vx) <= 2)
  expect_true(ncol(fit$vy) <= 2)

  # Project from X to latent
  Zx <- project(fit, X, source = "X")
  expect_equal(dim(Zx), c(n, ncol(fit$vx)))  # n x #components

  # Transfer from X to Y domain
  Yhat <- transfer(fit, X, source = "X", target = "Y")
  expect_equal(dim(Yhat), c(n, q))
})

test_that("rpls with 'ridge' penalty returns correct class and dims", {
  set.seed(123)
  n <- 15
  p <- 5
  q <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)
  X <- scale(X, center = TRUE, scale = FALSE)
  Y <- scale(Y, center = TRUE, scale = FALSE)

  # Fit a ridge RPLS
  fit_ridge <- rpls(X, Y, K = 2, lambda = 0.2, penalty = "ridge")

  # Basic checks
  expect_s3_class(fit_ridge, c("rpls", "cross_projector", "projector"))
  expect_true(ncol(fit_ridge$vx) <= 2)
  expect_true(ncol(fit_ridge$vy) <= 2)

  # Check projection from Y
  Zy <- project(fit_ridge, Y, source = "Y")
  expect_equal(dim(Zy), c(n, ncol(fit_ridge$vy)))

  # Transfer from Y -> X
  Xhat <- multivarious:::transfer(fit_ridge, Y, source = "Y", target = "X")
  expect_equal(dim(Xhat), c(n, p))
})


## tests/testthat/test-rpls.R
## (This can be appended to your existing test file or placed in a new file.)

test_that("rpls with user-supplied Q (GPLS) works for a small example", {
  skip_if_not_installed("Matrix")  # if using a sparse/dense Q
  library(Matrix)

  set.seed(999)
  n <- 12
  p <- 5
  q <- 3

  # Generate X, Y
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)
  # Center them
  X <- scale(X, center = TRUE, scale = FALSE)
  Y <- scale(Y, center = TRUE, scale = FALSE)

  # Create a random PSD matrix Q by crossproduct:
  Q0 <- matrix(rnorm(p * p), p, p)
  Qmat <- crossprod(Q0)  # p x p
  # Optionally, scale it so it doesn't blow up norms:
  Qmat <- Qmat / (1e3 + norm(Qmat, "F"))

  # Fit RPLS with Q => GPLS. Use ridge penalty as a test
  fit_gpls <- rpls(X, Y, K = 2, lambda = 0.05, penalty = "ridge",
                   Q = Qmat, verbose = FALSE)

  expect_s3_class(fit_gpls, c("rpls", "cross_projector", "projector"))
  # Check dimensions of loadings
  expect_true(ncol(fit_gpls$vx) <= 2)
  expect_true(ncol(fit_gpls$vy) <= 2)

  # Project X -> latent
  Zx <- project(fit_gpls, X, source = "X")
  expect_equal(dim(Zx), c(n, ncol(fit_gpls$vx)))

  # Transfer X -> Y
  Yhat <- transfer(fit_gpls, X, source = "X", target = "Y")
  expect_equal(dim(Yhat), c(n, q))

  # Just check we didn't blow up numerically
  expect_false(anyNA(Yhat))
})

test_that("rpls with nonnegative l1 penalty produces nonnegative loadings", {
  set.seed(234)
  n <- 10
  p <- 6
  q <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)

  X <- scale(X, center = TRUE, scale = FALSE)
  Y <- scale(Y, center = TRUE, scale = FALSE)

  # Nonnegative lasso
  fit_nonneg <- rpls(X, Y, K = 2, lambda = 0.1,
                     penalty = "l1", nonneg = TRUE,
                     verbose = FALSE)
  expect_s3_class(fit_nonneg, c("rpls", "cross_projector", "projector"))

  # The X-loadings (vx) should all be >= 0
  vx_mat <- fit_nonneg$vx
  expect_true(all(vx_mat >= 0))

  # Also check that we can project and transfer
  Fx <- project(fit_nonneg, X, source = "X")
  expect_equal(dim(Fx), c(n, ncol(vx_mat)))

  Yhat <- transfer(fit_nonneg, X, source = "X", target = "Y")
  expect_equal(dim(Yhat), c(n, q))
})

test_that("rpls partial projection works if partial cols are requested", {
  set.seed(567)
  n <- 12
  p <- 6
  q <- 4
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)
  X <- scale(X, center = TRUE, scale = FALSE)
  Y <- scale(Y, center = TRUE, scale = FALSE)

  # Fit standard l1-lasso RPLS
  fit_l1 <- rpls(X, Y, K = 3, lambda = 0.2, penalty = "l1", verbose = FALSE)

  # Suppose we only want partial columns 1:2
  part_cols <- 1:2

  # Subset X to these columns:
  X_sub <- X[, part_cols, drop = FALSE]   # Now X_sub is (12 x 2)

  # partial_project from X using X_sub and the same part_cols
  skip_if_not("partial_project.cross_projector" %in%
                as.character(methods(partial_project)))

  Fx_part <- partial_project(fit_l1,
                             new_data = X_sub,
                             colind   = part_cols,
                             source   = "X")
  # partial_project returns projection to full K-dimensional space

  expect_equal(dim(Fx_part), c(n, fit_l1$ncomp))

  # (You can do further checks, e.g. partial inverse projection, etc.)
})
