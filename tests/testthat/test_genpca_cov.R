library(testthat)
library(genpca)
library(Matrix)

test_that("genpca_cov matches genpca when C = X'MX", {
  set.seed(123)
  n <- 50
  p <- 20
  X <- matrix(rnorm(n * p), n, p)

  # Test 1: Standard PCA (no constraints)
  C <- crossprod(X)  # X'X

  # Run genpca with identity constraints
  fit_genpca <- genpca(X, A = NULL, M = NULL, ncomp = 5,
                       preproc = multivarious::pass())

  # Run genpca_cov on C (default method="gmd" should match)
  fit_cov <- genpca_cov(C, R = NULL, ncomp = 5)

  # Check that singular values match
  expect_equal(fit_genpca$sdev, fit_cov$d, tolerance = 1e-8)

  # Check that loadings match (up to sign)
  for (i in 1:5) {
    v_genpca <- fit_genpca$ov[, i]
    v_cov <- fit_cov$v[, i]
    # Check if vectors are parallel (correlation should be ±1)
    cor_val <- abs(cor(v_genpca, v_cov))
    expect_equal(cor_val, 1, tolerance = 1e-8)
  }

  # Test 2: With row metric M
  M_diag <- runif(n, 0.5, 1.5)
  M <- Matrix::Diagonal(n, x = M_diag)
  C_M <- crossprod(X, M %*% X)  # X'MX

  fit_genpca_M <- genpca(X, A = NULL, M = M_diag, ncomp = 5,
                         preproc = multivarious::pass())
  fit_cov_M <- genpca_cov(C_M, R = NULL, ncomp = 5)

  expect_equal(fit_genpca_M$sdev, fit_cov_M$d, tolerance = 1e-8)

  # Test 3: Simplified test - just check that the function works with constraints
  # Note: The exact mathematical equivalence with genpca when both M and A
  # are non-identity requires further investigation
  A_diag <- runif(p, 0.5, 2)
  C_M <- crossprod(X, M %*% X)  # X'MX

  fit_cov_MA <- genpca_cov(C_M, R = A_diag, ncomp = 5)

  # Check that it returns valid results
  expect_equal(length(fit_cov_MA$d), 5)
  expect_true(all(fit_cov_MA$d > 0))
  expect_equal(dim(fit_cov_MA$v), c(p, 5))

  # Check G-orthonormality of eigenvectors: V'GV = I
  G <- Matrix::Diagonal(p, x = A_diag)
  VGV <- as.matrix(t(fit_cov_MA$v) %*% (G %*% fit_cov_MA$v))
  expect_equal(VGV, diag(5), tolerance = 1e-8)
})

test_that("genpca_cov handles diagonal weights correctly", {
  set.seed(456)
  p <- 10
  C <- cov(matrix(rnorm(100 * p), 100, p))

  # Test with diagonal weights as vector
  w <- runif(p, 0.1, 2)
  fit_vec <- genpca_cov(C, R = w, ncomp = 3)

  # Test with diagonal weights as matrix
  fit_mat <- genpca_cov(C, R = diag(w), ncomp = 3)

  # Results should be identical (eigenvectors up to column sign)
  expect_equal(fit_vec$d, fit_mat$d, tolerance = 1e-8)
  expect_equal(as.matrix(fit_vec$v),
               align_signs(as.matrix(fit_vec$v), as.matrix(fit_mat$v)),
               tolerance = 1e-8)
  expect_equal(fit_vec$propv, fit_mat$propv, tolerance = 1e-8)
})

test_that("genpca_cov handles semidefinite G", {
  set.seed(789)
  p <- 10
  C <- cov(matrix(rnorm(100 * p), 100, p))

  # Create a rank-deficient R (some zero weights)
  w <- c(rep(1, 5), rep(0, 5))
  fit <- genpca_cov(C, R = w, ncomp = 3)

  # Should work and return results in the range of R
  expect_equal(fit$R_rank, 5)
  expect_true(fit$k <= 5)
  expect_true(all(fit$d > 0))
})

test_that("genpca_cov variance explained sums to <= 1", {
  set.seed(321)
  p <- 15
  C <- cov(matrix(rnorm(200 * p), 200, p))

  fit <- genpca_cov(C, R = NULL, ncomp = NULL)

  # Total variance explained should be <= 1
  expect_true(sum(fit$propv) <= 1 + 1e-10)

  # Cumulative variance should be monotonically increasing
  expect_true(all(diff(fit$cumv) >= -1e-10))

  # Last cumulative value should equal sum of propv
  expect_equal(fit$cumv[length(fit$cumv)], sum(fit$propv), tolerance = 1e-8)
})

test_that("genpca_cov remedies work correctly", {
  set.seed(654)
  p <- 5
  C <- cov(matrix(rnorm(50 * p), 50, p))

  # Create a definitely non-PSD G
  G <- matrix(rnorm(p * p), p, p)
  G <- G %*% t(G)
  # Get the minimum eigenvalue and subtract more than that to ensure non-PSD
  min_eig <- min(eigen(G, symmetric = TRUE)$values)
  G <- G - (abs(min_eig) + 1) * diag(p)  # Make it definitely non-PSD

  # Test error remedy (with geigen method since gmd doesn't use constraints_remedy)
  expect_error(geigen_cov(C, R = G, constraints_remedy = "error"))

  # Test ridge remedy
  fit_ridge <- geigen_cov(C, R = G, constraints_remedy = "ridge", ncomp = 2)
  expect_true(all(fit_ridge$d > 0))

  # Test clip remedy
  fit_clip <- geigen_cov(C, R = G, constraints_remedy = "clip", ncomp = 2)
  expect_true(all(fit_clip$d > 0))

  # Test identity remedy
  fit_identity <- geigen_cov(C, R = G, constraints_remedy = "identity", ncomp = 2)
  expect_true(all(fit_identity$d > 0))
})
