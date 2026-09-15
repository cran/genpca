testthat::test_that("genpls matches gplssvd_op on centered identity metrics", {
  skip_if_not_installed("Matrix")
  library(Matrix)
  set.seed(7)
  N <- 40
  I <- 12
  J <- 10
  k <- 3
  X <- matrix(rnorm(N * I), N, I)
  Y <- matrix(rnorm(N * J), N, J)

  # Operator path (centered)
  op <- gplssvd_op(X, Y, XLW = NULL, YLW = NULL, XRW = NULL, YRW = NULL,
                   k = k, center = TRUE, scale = FALSE)

  # High-level path with equivalent preprocessing
  fit <- genpls(X, Y, ncomp = k,
                preproc_x = multivarious::center(),
                preproc_y = multivarious::center())

  # Helpers
  align_signs <- function(A_ref, A_test) {
    A_ref <- as.matrix(A_ref)
    A_test <- as.matrix(A_test)
    ss <- sign(colSums(A_ref * A_test))
    ss[is.na(ss) | ss == 0] <- 1
    A_test %*% diag(ss, ncol(A_test))
  }
  expect_close <- function(A, B, tol = 1e-6) {
    A <- as.matrix(A)
    B <- as.matrix(B)
    testthat::expect_lt(norm(A - B, "F"), tol * (1 + norm(A, "F")))
  }

  # Compare spectra
  testthat::expect_equal(fit$d, op$d, tolerance = 1e-6)

  # p, q, fi, fj, lx, ly up to sign
  expect_close(align_signs(op$p,  fit$p),  op$p,  tol = 1e-6)
  expect_close(align_signs(op$q,  fit$q),  op$q,  tol = 1e-6)
  expect_close(align_signs(op$fi, fit$fi), op$fi, tol = 1e-6)
  expect_close(align_signs(op$fj, fit$fj), op$fj, tol = 1e-6)
  expect_close(align_signs(op$lx, fit$lx), op$lx, tol = 1e-6)
  expect_close(align_signs(op$ly, fit$ly), op$ly, tol = 1e-6)
})

testthat::test_that("genpls matches gplssvd_op with diagonal row metrics", {
  skip_if_not_installed("Matrix")
  library(Matrix)
  set.seed(8)
  N <- 40
  I <- 10
  J <- 9
  k <- 3
  X <- matrix(rnorm(N * I), N, I)
  Y <- matrix(rnorm(N * J), N, J)
  wX <- runif(N)
  wX <- wX / sum(wX)
  wY <- runif(N)
  wY <- wY / sum(wY)
  MX <- Diagonal(x = wX)
  MY <- Diagonal(x = wY)

  op <- gplssvd_op(X, Y, XLW = MX, YLW = MY, XRW = NULL, YRW = NULL,
                   k = k, center = TRUE, scale = FALSE)

  fit <- genpls(X, Y, Mx = MX, My = MY, ncomp = k,
                preproc_x = multivarious::center(),
                preproc_y = multivarious::center())

  align_signs <- function(A_ref, A_test) {
    A_ref <- as.matrix(A_ref)
    A_test <- as.matrix(A_test)
    ss <- sign(colSums(A_ref * A_test))
    ss[is.na(ss) | ss == 0] <- 1
    A_test %*% diag(ss, ncol(A_test))
  }
  expect_close <- function(A, B, tol = 1e-6) {
    A <- as.matrix(A)
    B <- as.matrix(B)
    testthat::expect_lt(norm(A - B, "F"), tol * (1 + norm(A, "F")))
  }

  testthat::expect_equal(fit$d, op$d, tolerance = 1e-6)
  expect_close(align_signs(op$p,  fit$p),  op$p,  tol = 1e-6)
  expect_close(align_signs(op$q,  fit$q),  op$q,  tol = 1e-6)
  expect_close(align_signs(op$fi, fit$fi), op$fi, tol = 1e-6)
  expect_close(align_signs(op$fj, fit$fj), op$fj, tol = 1e-6)
  expect_close(align_signs(op$lx, fit$lx), op$lx, tol = 1e-6)
  expect_close(align_signs(op$ly, fit$ly), op$ly, tol = 1e-6)
})
