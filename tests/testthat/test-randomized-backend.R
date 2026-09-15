test_that("randomized backend returns metric-orthonormal factors", {
  set.seed(101)
  n <- 120
  p <- 180
  k <- 8

  X <- matrix(rnorm(n * p), n, p)
  Q <- Matrix::Diagonal(n, x = runif(n, 0.8, 1.3))
  R <- Matrix::Diagonal(p, x = runif(p, 0.7, 1.4))

  fit <- genpca(
    X,
    M = Q,
    A = R,
    ncomp = k,
    method = "randomized",
    oversample = 15L,
    n_power = 1L,
    n_polish = 1L,
    preproc = multivarious::pass()
  )

  QtU <- crossprod(fit$ou, Q %*% fit$ou)
  RtV <- crossprod(fit$ov, R %*% fit$ov)

  expect_lt(max(abs(QtU - diag(ncol(QtU)))), 1e-6)
  expect_lt(max(abs(RtV - diag(ncol(RtV)))), 1e-6)
})

test_that("randomized backend tracks leading eigen variance on moderate problems", {
  set.seed(202)
  n <- 140
  p <- 100
  k <- 10

  X <- matrix(rnorm(n * p), n, p)
  Q <- Matrix::Diagonal(n, x = runif(n, 0.9, 1.2))
  R <- Matrix::Diagonal(p, x = runif(p, 0.9, 1.2))

  fit_eig <- genpca(
    X,
    M = Q,
    A = R,
    ncomp = k,
    method = "eigen",
    preproc = multivarious::pass()
  )

  fit_rand <- genpca(
    X,
    M = Q,
    A = R,
    ncomp = k,
    method = "randomized",
    oversample = 20L,
    n_power = 1L,
    n_polish = 1L,
    preproc = multivarious::pass()
  )

  var_ratio <- sum(fit_rand$sdev^2) / sum(fit_eig$sdev^2)

  expect_true(is.finite(var_ratio))
  expect_gt(var_ratio, 0.9)
  expect_equal(fit_rand$method, "randomized")
})

test_that("randomized polish tolerance validates and runs", {
  set.seed(303)
  n <- 100
  p <- 140
  k <- 6

  X <- matrix(rnorm(n * p), n, p)
  Q <- Matrix::Diagonal(n, x = runif(n, 0.8, 1.4))
  R <- Matrix::Diagonal(p, x = runif(p, 0.8, 1.4))

  expect_error(
    genpca(
      X,
      M = Q,
      A = R,
      ncomp = k,
      method = "randomized",
      n_polish = 2L,
      tol_polish_randomized = -1,
      preproc = multivarious::pass()
    ),
    "tol_polish_randomized must be a single non-negative number."
  )

  fit <- genpca(
    X,
    M = Q,
    A = R,
    ncomp = k,
    method = "randomized",
    oversample = 10L,
    n_power = 1L,
    n_polish = 3L,
    tol_polish_randomized = 1e-4,
    preproc = multivarious::pass()
  )

  expect_equal(fit$method, "randomized")
  expect_true(all(is.finite(fit$sdev)))
  expect_equal(length(fit$sdev), k)
})
