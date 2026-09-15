test_that("randomized full sketches preserve the analytic spectrum across seeds", {
  # Sign sketches collide for many seeds; std::uniform_int_distribution maps
  # engines differently under libc++ and libstdc++, so one fixed seed hid it.
  X <- diag(c(1, 1e6)); W <- diag(c(1, 1e-10))
  for (cpp in c(FALSE, TRUE)) for (seed in c(1L, 2L, 1234L)) {
    for (side in c("rows", "columns")) for (jitter in c(0, 1e-10)) {
      M <- if (side == "rows") W else diag(2)
      A <- if (side == "columns") W else diag(2)
      f <- genpca:::gmd_randomized(X, M, A, k = 2, seed = seed,
                                   jitter = jitter, use_cpp = cpp, n_polish = 1)
      expect_equal(f$d, c(10, 1), tolerance = 1e-7)
      expect_equal(crossprod(f$u, M %*% f$u), diag(2), tolerance = 1e-7)
      expect_equal(crossprod(f$v, A %*% f$v), diag(2), tolerance = 1e-7)
      expect_equal(sweep(f$u, 2, f$d, "*") %*% t(f$v), X, tolerance = 1e-7)
    }
  }
})

test_that("Gram jitter never fabricates independent metric directions", {
  # A rank-one sketch must yield one genuinely normalized direction.
  for (jitter in c(0, 1e-10, 1e-3)) {
    B <- genpca:::metric_orthonormalize(cbind(c(1, 1e6), c(1, 1e6)),
                                        function(x) x, jitter = jitter)
    expect_equal(ncol(B), 1L)
    expect_equal(crossprod(B), matrix(1), tolerance = 1e-10)
  }
  # Exercise the corresponding C++ path through a rank-deficient data block,
  # with a full column sketch but fewer columns than rows.
  X <- cbind(c(1, 1e6, 0), c(1, 1e6, 0))
  for (cpp in c(FALSE, TRUE)) {
    f <- genpca:::gmd_randomized(X, diag(3), diag(2), 2, use_cpp = cpp)
    expect_length(f$d, 1L)
    expect_equal(f$d, sqrt(2 * (1 + 1e12)), tolerance = 1e-10)
    expect_equal(crossprod(f$u), matrix(1), tolerance = 1e-10)
    expect_equal(as.matrix(crossprod(f$v)), matrix(1), tolerance = 1e-10)
  }
})

test_that("MLE rescale delta measures only the reciprocal scaling penalty", {
  set.seed(9)
  X <- matrix(rnorm(60), 10, 6)
  lambda <- 0.02
  raw <- gpca_mle(X, ncomp = 2, max_iter = 3, lambda = lambda, scale_fix = "none")
  penalty <- function(f) lambda * (ncol(X) * sum(diag(f$M)) +
                                   nrow(X) * sum(diag(f$A)))
  expect_identical(raw$loglik_rescale_delta, 0)
  for (sf in c("trace", "det")) {
    f <- gpca_mle(X, ncomp = 2, max_iter = 3, lambda = lambda, scale_fix = sf)
    expect_equal(f$loglik_rescale_delta, -0.5 * (penalty(f) - penalty(raw)),
                 tolerance = 1e-10)
    expect_equal(f$loglik, tail(f$loglik_path, 1) +
                   f$loglik_rescale_delta + f$loglik_refit_delta, tolerance = 1e-12)
    # Independent direct precision determinant formulation of the objective.
    E <- X - multivarious::reconstruct(f$fit)
    logdet <- function(M) as.numeric(determinant(as.matrix(M), logarithm = TRUE)$modulus)
    ll <- 0.5 * (ncol(X) * logdet(f$M) + nrow(X) * logdet(f$A) -
                  sum(E * as.matrix(f$M %*% E %*% f$A)) - penalty(f))
    expect_equal(f$loglik, ll, tolerance = 1e-8)
  }
})

test_that("partial Gaussian sketches recover a known weighted low-rank matrix", {
  set.seed(14)
  n <- 12; p <- 8
  U <- qr.Q(qr(matrix(rnorm(n * 3), n)))
  V <- qr.Q(qr(matrix(rnorm(p * 3), p)))
  q <- 10^seq(-5, 0, length.out = n)
  r <- 10^seq(-5, 0, length.out = p)
  # Analytic whitened SVD, independent of the solver's Gram calculations.
  X <- sweep(U, 1, sqrt(q), "/") %*% diag(c(4, 2, .5)) %*%
    t(sweep(V, 1, sqrt(r), "/"))
  for (cpp in c(FALSE, TRUE)) for (seed in c(1L, 1234L)) {
    for (sparse_q in c(FALSE, TRUE)) for (sparse_r in c(FALSE, TRUE)) {
      Q <- if (sparse_q) Matrix::Diagonal(x = q) else diag(q)
      R <- if (sparse_r) Matrix::Diagonal(x = r) else diag(r)
      f <- genpca:::gmd_randomized(X, Q, R, 3, oversample = 2,
                                   n_power = 1, n_polish = 1, seed = seed, use_cpp = cpp)
      expect_equal(f$d, c(4, 2, .5), tolerance = 1e-7)
      expect_equal(as.matrix(crossprod(f$u, Q %*% f$u)), diag(3), tolerance = 1e-7)
      expect_equal(as.matrix(crossprod(f$v, R %*% f$v)), diag(3), tolerance = 1e-7)
      E <- X - sweep(f$u, 2, f$d, "*") %*% t(f$v)
      expect_lt(sqrt(sum(E * as.matrix(Q %*% E %*% R))), 1e-7)
    }
  }
})
