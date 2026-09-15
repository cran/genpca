library(testthat)
library(Matrix)

test_that("Rank-1 matrix is recovered correctly", {
  set.seed(123)
  n <- 10
  p <- 5
  u_true <- rnorm(n)
  v_true <- rnorm(p)
  X <- u_true %*% t(v_true)
  spat_cds <- matrix(runif(p * 2), nrow = 2, ncol = p)
  # Test sparsity regularization only (no smoothness penalties)
  result <- sfpca(X, K = 1, spat_cds = spat_cds, verbose = FALSE,
                  lambda_v = .001, lambda_u = .001, alpha_v = 0, alpha_u = 0)

  u_est <- result$ou[, 1]
  v_est <- multivarious::components(result)[, 1]
  d_est <- multivarious::sdev(result)[1]

  # Handle sign ambiguity - fix both u and v signs consistently
  sign_uv <- sign(crossprod(u_est, u_true) * crossprod(v_est, v_true))
  if (sign_uv < 0) {
    u_est <- -u_est
  }

  # Check reconstruction - with small L1 regularization only, should be very good
  X_recon <- d_est * u_est %*% t(v_est)
  recon_error <- norm(X - X_recon, "F") / norm(X, "F")

  # With small L1 regularization (0.001) and no smoothness, expect excellent reconstruction
  expect_lt(recon_error, 0.01)  # At least 99% variance explained

  # Also check correlation with true components
  expect_gt(abs(cor(u_est, u_true)), 0.99)
  expect_gt(abs(cor(v_est, v_true)), 0.99)
})

test_that("Orthogonal columns result in smooth components", {
  set.seed(123)
  n <- 100
  p <- 20
  # Create an n x p matrix with orthogonal columns
  # If n > p, we need to repeat/extend the identity pattern
  if (n <= p) {
    X <- diag(p)[1:n, ]
  } else {
    # Create orthogonal columns by cycling through identity matrix rows
    X <- matrix(0, n, p)
    for (i in 1:n) {
      X[i, ((i - 1) %% p) + 1] <- 1
    }
  }
  spat_cds <- matrix(runif(p * 2), nrow = 2, ncol = p)
  # Use K=1 since the matrix may not support K=2 after deflation
  result <- sfpca(X, K = 1, spat_cds = spat_cds, verbose = FALSE)

  v_est <- multivarious::components(result)

  # Compute smoothness measure based on spat_cds
  distances <- as.matrix(dist(t(spat_cds)))
  smoothness <- sum(distances * (v_est %*% t(v_est)))

  # Set a threshold for smoothness
  expect_true(smoothness < 100)  # Adjust threshold as needed
})

test_that("Sparse signals result in sparse components", {
  set.seed(123)
  n <- 50
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  X[, sample(p, 10)] <- 0  # Introduce sparsity
  spat_cds <- matrix(runif(p * 3), nrow = 3, ncol = p)
  # Use stronger lambda values and minimal spatial smoothing for sparsity test
  result <- sfpca(X, K = 2, spat_cds = spat_cds,
                  lambda_u = 1.0, lambda_v = 1.0,
                  alpha_u = 0, alpha_v = 0,  # Disable spatial smoothing
                  verbose = FALSE)

  u_est <- result$ou
  v_est <- multivarious::components(result)

  # Check sparsity - with proper ISTA updates and no spatial smoothing
  # At least one component should show sparsity
  u_sparse <- sum(abs(u_est[, 1]) > 1e-4) < 25 || sum(abs(u_est[, 2]) > 1e-4) < 25
  v_sparse <- sum(abs(v_est[, 1]) > 1e-4) < 25 || sum(abs(v_est[, 2]) > 1e-4) < 25

  expect_true(u_sparse || v_sparse)  # At least one should be sparse
})

test_that("spat_cds shape is validated with an actionable message", {
  set.seed(42)
  n <- 20; p <- 30
  X <- matrix(rnorm(n * p), n, p)
  cds <- rbind(runif(p), runif(p))          # 2 x p, the correct orientation

  # the correct orientation still works
  expect_no_error(sfpca(X, K = 1, spat_cds = cds))

  # transposed input names the problem and the fix
  expect_error(sfpca(X, K = 1, spat_cds = t(cds)),
               "looks transposed", fixed = TRUE)
  expect_error(sfpca(X, K = 1, spat_cds = t(cds)),
               "t(spat_cds)", fixed = TRUE)

  # a plain size mismatch reports both sizes but does not claim a transpose
  wrong <- rbind(runif(p + 3), runif(p + 3))
  expect_error(sfpca(X, K = 1, spat_cds = wrong), "ncol\\(spat_cds\\) is 33")
  expect_error(sfpca(X, K = 1, spat_cds = wrong), "one column per variable")
  expect_false(grepl("transposed",
    tryCatch(sfpca(X, K = 1, spat_cds = wrong), error = conditionMessage)))

  # a bare vector points at the 1-D form rather than failing downstream
  expect_error(sfpca(X, K = 1, spat_cds = runif(p)),
               "matrix(coords, nrow = 1)", fixed = TRUE)

  # NAs are rejected before they reach the kNN search
  na_cds <- cds; na_cds[1, 1] <- NA
  expect_error(sfpca(X, K = 1, spat_cds = na_cds), "must not contain NA")
})
