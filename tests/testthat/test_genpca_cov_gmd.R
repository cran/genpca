library(testthat)
library(genpca)
library(Matrix)

test_that("genpca_cov with method='gmd' exactly matches two-sided genpca", {
  set.seed(123)
  n <- 50
  p <- 20
  k <- 5
  X <- matrix(rnorm(n * p), n, p)

  # Test 1: Identity constraints (standard PCA)
  C <- crossprod(X)  # X'X

  fit_genpca <- genpca(X, A = NULL, M = NULL, ncomp = k,
                       preproc = multivarious::pass())
  fit_cov_gmd <- genpca_cov(C, R = NULL, ncomp = k)

  # Should match exactly
  expect_equal(fit_genpca$sdev, fit_cov_gmd$d, tolerance = 1e-8,
               label = "GMD: Singular values with identity constraints")

  # Eigenvectors should match (up to sign)
  for (i in 1:k) {
    cor_val <- abs(cor(fit_genpca$ov[, i], fit_cov_gmd$v[, i]))
    expect_equal(cor_val, 1, tolerance = 1e-8,
                 label = paste("GMD: Eigenvector", i, "correlation (identity)"))
  }
})

test_that("genpca_cov GMD matches genpca with column constraint A", {
  set.seed(456)
  n <- 40
  p <- 15
  k <- 5
  X <- matrix(rnorm(n * p), n, p)

  # Column constraint A (this is where geigen failed)
  A_diag <- runif(p, 0.5, 2)
  C <- crossprod(X)  # X'X (M = I)

  fit_genpca <- genpca(X, A = A_diag, M = NULL, ncomp = k,
                       preproc = multivarious::pass())
  fit_cov_gmd <- genpca_cov(C, R = A_diag, ncomp = k)

  # Now with GMD method, these should match exactly!
  expect_equal(fit_genpca$sdev, fit_cov_gmd$d, tolerance = 1e-8,
               label = "GMD: Singular values with A constraint")

  # Check eigenvectors match (up to sign)
  for (i in 1:k) {
    cor_val <- abs(cor(fit_genpca$ov[, i], fit_cov_gmd$v[, i]))
    expect_equal(cor_val, 1, tolerance = 1e-8,
                 label = paste("GMD: Eigenvector", i, "with A constraint"))
  }

  # Verify A-orthonormality
  A_mat <- diag(A_diag)
  VAV <- t(fit_cov_gmd$v) %*% A_mat %*% fit_cov_gmd$v
  expect_equal(as.matrix(VAV), diag(k), tolerance = 1e-8,
               label = "GMD: A-orthonormality of eigenvectors")
})

test_that("genpca_cov GMD matches genpca with row constraint M", {
  set.seed(789)
  n <- 30
  p <- 10
  k <- 5
  X <- matrix(rnorm(n * p), n, p)

  # Row constraint M, column constraint A = I
  M_diag <- runif(n, 0.5, 1.5)
  M <- Matrix::Diagonal(n, x = M_diag)
  C_M <- crossprod(X, M %*% X)  # X'MX

  fit_genpca <- genpca(X, A = NULL, M = M_diag, ncomp = k,
                       preproc = multivarious::pass())
  fit_cov_gmd <- genpca_cov(C_M, R = NULL, ncomp = k)

  # Should match
  expect_equal(fit_genpca$sdev, fit_cov_gmd$d, tolerance = 1e-8,
               label = "GMD: Singular values with M constraint")

  # Check eigenvectors
  for (i in 1:k) {
    cor_val <- abs(cor(fit_genpca$ov[, i], fit_cov_gmd$v[, i]))
    expect_equal(cor_val, 1, tolerance = 1e-8,
                 label = paste("GMD: Eigenvector", i, "with M constraint"))
  }
})

test_that("genpca_cov GMD matches genpca with both M and A constraints", {
  set.seed(321)
  n <- 40
  p <- 15
  k <- 5
  X <- matrix(rnorm(n * p), n, p)

  # Both row and column constraints
  M_diag <- runif(n, 0.5, 1.5)
  A_diag <- runif(p, 0.5, 2)

  M <- Matrix::Diagonal(n, x = M_diag)
  C_M <- crossprod(X, M %*% X)  # X'MX

  fit_genpca <- genpca(X, A = A_diag, M = M_diag, ncomp = k,
                       preproc = multivarious::pass())
  fit_cov_gmd <- genpca_cov(C_M, R = A_diag, ncomp = k)

  # Should match exactly with GMD method
  expect_equal(fit_genpca$sdev, fit_cov_gmd$d, tolerance = 1e-8,
               label = "GMD: Singular values with both M and A")

  # Check eigenvectors
  for (i in 1:k) {
    cor_val <- abs(cor(fit_genpca$ov[, i], fit_cov_gmd$v[, i]))
    expect_equal(cor_val, 1, tolerance = 1e-8,
                 label = paste("GMD: Eigenvector", i, "with both M and A"))
  }

  # Verify A-orthonormality
  A_mat <- diag(A_diag)
  VAV <- t(fit_cov_gmd$v) %*% A_mat %*% fit_cov_gmd$v
  expect_equal(as.matrix(VAV), diag(k), tolerance = 1e-8,
               label = "GMD: A-orthonormality with both constraints")
})

test_that("genpca_cov variance explained matches between GMD and two-sided", {
  set.seed(654)
  n <- 50
  p <- 20
  X <- matrix(rnorm(n * p), n, p)

  # Test with various constraint combinations
  test_cases <- list(
    list(M = NULL, A = NULL, desc = "Identity constraints"),
    list(M = runif(n, 0.5, 2), A = NULL, desc = "M constraint only"),
    list(M = NULL, A = runif(p, 0.5, 2), desc = "A constraint only"),
    list(M = runif(n, 0.5, 2), A = runif(p, 0.5, 2), desc = "Both constraints")
  )

  for (tc in test_cases) {
    # Compute covariance
    if (is.null(tc$M)) {
      C <- crossprod(X)
    } else {
      M <- Matrix::Diagonal(n, x = tc$M)
      C <- crossprod(X, M %*% X)
    }

    # Run both methods
    fit_genpca <- genpca(X, A = tc$A, M = tc$M, ncomp = 5,
                         preproc = multivarious::pass())
    fit_cov_gmd <- genpca_cov(C, R = tc$A, ncomp = 5)

    # Variance explained should match
    expect_equal(fit_genpca$propv, fit_cov_gmd$propv, tolerance = 1e-8,
                 label = paste("GMD propv:", tc$desc))
    expect_equal(fit_genpca$cumv, fit_cov_gmd$cumv, tolerance = 1e-8,
                 label = paste("GMD cumv:", tc$desc))
  }
})

test_that("genpca_cov correctly identifies method='geigen' vs 'gmd' differences", {
  set.seed(987)
  n <- 30
  p <- 15
  X <- matrix(rnorm(n * p), n, p)

  # With A constraint (where geigen and gmd differ)
  A_diag <- runif(p, 0.5, 2)
  C <- crossprod(X)

  fit_gmd <- genpca_cov(C, R = A_diag, ncomp = 5)
  fit_geigen <- geigen_cov(C, R = A_diag, ncomp = 5)

  # Methods should be recorded
  expect_equal(fit_gmd$method, "gmd")
  expect_equal(fit_geigen$method, "geigen")

  # Results should differ (as we discovered)
  expect_false(isTRUE(all.equal(fit_gmd$d, fit_geigen$d, tolerance = 1e-8)),
               label = "GMD and geigen produce different results with A constraint")

  # But both should produce A-orthonormal vectors
  A_mat <- diag(A_diag)

  VAV_gmd <- t(fit_gmd$v) %*% A_mat %*% fit_gmd$v
  expect_equal(as.matrix(VAV_gmd), diag(5), tolerance = 1e-8,
               label = "GMD produces A-orthonormal vectors")

  VAV_geigen <- t(fit_geigen$v) %*% A_mat %*% fit_geigen$v
  expect_equal(as.matrix(VAV_geigen), diag(5), tolerance = 1e-8,
               label = "geigen produces A-orthonormal vectors")
})
