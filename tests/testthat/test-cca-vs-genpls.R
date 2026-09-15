testthat::test_that("genpls reproduces canonical correlations on toy data (CCA via precisions)", {
  skip_if_not_installed("Matrix")
  library(Matrix)

  set.seed(42)
  n <- 150
  p <- 6
  q <- 4
  k <- min(p, q, 3)

  # Construct correlated blocks via shared latent factors
  r <- 3
  Z  <- matrix(rnorm(n * r), n, r)
  Bx <- matrix(rnorm(r * p), r, p)
  By <- matrix(rnorm(r * q), r, q)
  X  <- Z %*% Bx + 0.3 * matrix(rnorm(n * p), n, p)
  Y  <- Z %*% By + 0.3 * matrix(rnorm(n * q), n, q)

  # Center (no scaling) so covariances are meaningful
  Xc <- scale(X, center = TRUE, scale = FALSE)
  Yc <- scale(Y, center = TRUE, scale = FALSE)

  # Sample covariances and ridge-stabilized precisions
  Sx <- cov(Xc)
  Sy <- cov(Yc)
  lam <- 1e-8
  Ax <- solve(Sx + lam * diag(p))
  Ay <- solve(Sy + lam * diag(q))

  # CCA via stats::cancor
  cca <- stats::cancor(Xc, Yc, xcenter = FALSE, ycenter = FALSE)
  rho_cca <- cca$cor[seq_len(k)]

  # CCA via GPLSSVD mapping (pass centered data; column metrics are precisions)
  fit <- genplsc(Xc, Yc, Ax = Ax, Ay = Ay, ncomp = k,
                 preproc_x = multivarious::pass(),
                 preproc_y = multivarious::pass())

  # Canonical correlations are singular values of cross-correlation operator
  rho_g <- fit$d[seq_len(k)] / (n - 1)

  testthat::expect_equal(rho_g, rho_cca, tolerance = 1e-3)

  # Latent variable correlations should also match canonical correlations
  cc_lat <- sapply(seq_len(k), function(j) cor(fit$lx[, j], fit$ly[, j]))
  testthat::expect_equal(cc_lat, rho_cca, tolerance = 5e-2)
})
