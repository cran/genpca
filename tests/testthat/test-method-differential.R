test_that("spectra stays within a stable envelope versus eigen on diagonal metrics", {
  set.seed(812)

  # Case 1: n >= p
  n1 <- 140
  p1 <- 90
  k1 <- 12
  X1 <- matrix(rnorm(n1 * p1), n1, p1)
  M1 <- Matrix::Diagonal(n1, x = runif(n1, 0.8, 1.2))
  A1 <- Matrix::Diagonal(p1, x = runif(p1, 0.8, 1.2))

  fit_eig_1 <- genpca(
    X1, M = M1, A = A1, ncomp = k1,
    method = "eigen", preproc = multivarious::center()
  )
  fit_spc_1 <- genpca(
    X1, M = M1, A = A1, ncomp = k1,
    method = "spectra", preproc = multivarious::center()
  )

  rel1 <- max(abs(fit_spc_1$sdev - fit_eig_1$sdev) / pmax(abs(fit_eig_1$sdev), 1e-12))
  expect_lt(rel1, 1e-6)
  expect_lt(max(abs(crossprod(fit_spc_1$ou, M1 %*% fit_spc_1$ou) - diag(k1))), 1e-6)
  expect_lt(max(abs(crossprod(fit_spc_1$ov, A1 %*% fit_spc_1$ov) - diag(k1))), 1e-6)

  # Case 2: n < p
  n2 <- 90
  p2 <- 180
  k2 <- 10
  X2 <- matrix(rnorm(n2 * p2), n2, p2)
  M2 <- Matrix::Diagonal(n2, x = runif(n2, 0.8, 1.2))
  A2 <- Matrix::Diagonal(p2, x = runif(p2, 0.8, 1.2))

  fit_eig_2 <- genpca(
    X2, M = M2, A = A2, ncomp = k2,
    method = "eigen", preproc = multivarious::center()
  )
  fit_spc_2 <- genpca(
    X2, M = M2, A = A2, ncomp = k2,
    method = "spectra", preproc = multivarious::center()
  )

  rel2 <- max(abs(fit_spc_2$sdev - fit_eig_2$sdev) / pmax(abs(fit_eig_2$sdev), 1e-12))
  expect_lt(rel2, 1e-6)
  expect_lt(max(abs(crossprod(fit_spc_2$ou, M2 %*% fit_spc_2$ou) - diag(k2))), 1e-6)
  expect_lt(max(abs(crossprod(fit_spc_2$ov, A2 %*% fit_spc_2$ov) - diag(k2))), 1e-6)
})

test_that("randomized tracks eigen singular values within approximation envelope", {
  set.seed(813)
  n <- 180
  p <- 360
  k <- 20

  X <- matrix(rnorm(n * p), n, p)
  M <- Matrix::Diagonal(n, x = runif(n, 0.7, 1.3))
  A <- Matrix::Diagonal(p, x = runif(p, 0.7, 1.3))

  fit_eig <- genpca(
    X, M = M, A = A, ncomp = k,
    method = "eigen", preproc = multivarious::pass()
  )
  fit_rand <- genpca(
    X, M = M, A = A, ncomp = k,
    method = "randomized",
    oversample = 20L,
    n_power = 1L,
    n_polish = 1L,
    tol_polish_randomized = 1e-4,
    preproc = multivarious::pass()
  )

  rel <- max(abs(fit_rand$sdev - fit_eig$sdev) / pmax(abs(fit_eig$sdev), 1e-12))
  var_ratio <- sum(fit_rand$sdev^2) / sum(fit_eig$sdev^2)

  expect_lt(rel, 0.12)
  expect_gt(var_ratio, 0.9)
})

test_that("eigen backend factors sparse symmetric metrics exactly regardless of maxeig", {
  set.seed(814)
  n <- 140
  p <- 180
  k <- 12
  X <- matrix(rnorm(n * p), n, p)

  make_sparse_spd <- function(m) {
    methods::as(
      Matrix::Diagonal(m, x = rep(1.2, m)) +
        Matrix::bandSparse(
          m, m,
          k = c(-1, 1),
          diagonals = list(rep(-0.15, m - 1), rep(-0.15, m - 1))
        ),
      "dsCMatrix"
    )
  }

  M <- make_sparse_spd(n)
  A <- make_sparse_spd(p)

  fit <- genpca(
    X, M = M, A = A, ncomp = k,
    method = "eigen",
    maxeig = 80,
    preproc = multivarious::pass()
  )

  expect_equal(length(fit$sdev), k)
  expect_true(all(is.finite(fit$sdev)))
  fit_full <- genpca(
    X, M = M, A = A, ncomp = k,
    method = "eigen",
    maxeig = Inf,
    preproc = multivarious::pass()
  )
  expect_equal(fit$sdev, fit_full$sdev, tolerance = 1e-8)
})
