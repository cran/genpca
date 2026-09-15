test_that("deflation C++ handles dense metrics via sparse coercion", {
  set.seed(811)
  n <- 80
  p <- 50
  k <- 8

  X <- matrix(rnorm(n * p), n, p)
  M <- crossprod(matrix(rnorm(n * n), n, n)) / n + diag(n) * 0.5
  A <- crossprod(matrix(rnorm(p * p), p, p)) / p + diag(p) * 0.5

  fit_cpp <- genpca(
    X,
    M = M,
    A = A,
    ncomp = k,
    method = "deflation",
    use_cpp = TRUE,
    preproc = multivarious::center(),
    threshold = 1e-7,
    maxit_deflation = 400L
  )

  fit_r <- genpca(
    X,
    M = M,
    A = A,
    ncomp = k,
    method = "deflation",
    use_cpp = FALSE,
    preproc = multivarious::center(),
    threshold = 1e-7,
    maxit_deflation = 400L
  )

  expect_equal(fit_cpp$method, "deflation")
  expect_equal(length(fit_cpp$sdev), k)
  expect_true(all(is.finite(fit_cpp$sdev)))

  rel <- max(abs(fit_cpp$sdev - fit_r$sdev) / pmax(abs(fit_r$sdev), 1e-12))
  expect_lt(rel, 5e-2)
})
