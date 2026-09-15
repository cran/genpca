test_that("ensure_spd makes adjacency usable as a metric", {
  set.seed(1)
  n <- 200L
  idx <- sample(n * n, size = 1000L)
  A <- matrix(0, n, n)
  A[idx] <- 1
  A <- (A + t(A)) / 2
  diag(A) <- 0
  As <- Matrix::Matrix(A, sparse = TRUE)

  M <- ensure_spd(As)
  expect_true(inherits(M, "Matrix"))
  expect_silent(Matrix::Cholesky(M, LDL = FALSE, super = TRUE))
})

test_that("ensure_spd preserves already SPD matrices", {
  set.seed(42)
  n <- 50
  # Create a definitely SPD matrix
  B <- matrix(rnorm(n * n), n, n)
  S <- crossprod(B) + diag(n)  # SPD by construction

  S_spd <- ensure_spd(S)
  expect_true(inherits(S_spd, "Matrix"))
  # Should be minimal changes
  expect_equal(as.matrix(S_spd), S, tolerance = 1e-10)
})

test_that("ensure_spd handles small negative eigenvalues", {
  set.seed(123)
  n <- 30
  # Create a nearly SPD matrix with small negative eigenvalue
  B <- matrix(rnorm(n * n), n, n)
  S <- (B + t(B)) / 2
  E <- eigen(S)
  # Make smallest eigenvalue slightly negative
  E$values[n] <- -0.01
  S_bad <- E$vectors %*% diag(E$values) %*% t(E$vectors)

  S_fixed <- ensure_spd(S_bad)
  expect_true(inherits(S_fixed, "Matrix"))

  # Check it's now SPD via successful Cholesky
  expect_silent(Matrix::Cholesky(S_fixed, LDL = FALSE, super = TRUE))
})
