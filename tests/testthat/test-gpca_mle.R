library(testthat)
library(genpca)

test_that("gpca_mle returns SPD metrics with correct dimensions", {
  set.seed(1)
  X <- matrix(rnorm(60), 10, 6)
  res <- suppressWarnings(
    gpca_mle(X, ncomp = 3, max_iter = 4, lambda = 1e-3, verbose = FALSE)
  )
  expect_equal(dim(res$A), c(6, 6))
  expect_equal(dim(res$M), c(10, 10))
  expect_true(genpca:::is_spd(res$A))
  expect_true(genpca:::is_spd(res$M))
  expect_s3_class(res$fit, "genpca")
  expect_true(length(res$loglik_path) >= 1)
  expect_true(is.finite(res$loglik))
})

test_that("gpca_mle produces a finite likelihood path", {
  set.seed(2)
  X <- matrix(rnorm(40), 8, 5)
  res <- suppressWarnings(
    gpca_mle(X, ncomp = 2, max_iter = 5, lambda = 1e-3, verbose = FALSE)
  )
  expect_gt(length(res$loglik_path), 1)
  expect_true(all(is.finite(res$loglik_path)))
})

test_that("gpca_mle path stays monotone at large covariance scale (no tolerant ridge on Sigma)", {
  # Regression: the covariance updates Sigma_r = E A E'/p + lambda I have
  # minimum eigenvalue lambda; a relative PD margin used to shift them by a
  # multiple of their largest diagonal once max(diag) > 1e6 * lambda.
  set.seed(2)
  X <- matrix(rnorm(180), 30, 6)
  for (Xs in list(X, 100 * X)) {
    res <- gpca_mle(Xs, ncomp = 2, max_iter = 8)
    path <- res$loglik_path
    slack <- 1e-8 * pmax(abs(path[-length(path)]), 1)
    expect_true(all(diff(path) >= -slack))
    expect_true(genpca:::is_pd(res$M, rtol = 0))
    expect_true(genpca:::is_pd(res$A, rtol = 0))
    # No rescaling has exactly zero effect. The final refit/reevaluation
    # discrepancy is retained separately, not clamped away.
    expect_identical(res$loglik_rescale_delta, 0)
    expect_equal(res$loglik_refit_delta, res$loglik - tail(path, 1))
    expect_true(is.finite(res$loglik_refit_delta))
  }
})
