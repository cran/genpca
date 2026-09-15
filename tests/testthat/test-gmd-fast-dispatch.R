test_that("gmd_fast_cpp dense/sparse dispatch works", {
  set.seed(42)
  n <- 100
  p <- 30
  k <- 5
  X <- matrix(rnorm(n * p), n, p)
  Qd <- diag(runif(n, 0.5, 2))
  Rd <- diag(runif(p, 0.5, 2))

  Qs <- methods::as(Matrix::Matrix(Qd, sparse = TRUE), "dgCMatrix")
  Rs <- methods::as(Matrix::Matrix(Rd, sparse = TRUE), "dgCMatrix")

  # Test dense version
  res_dense <- genpca:::gmd_fast_cpp(X, Qd, Rd, k)
  expect_equal(length(res_dense$d), k)
  expect_equal(nrow(res_dense$u), n)
  expect_equal(nrow(res_dense$v), p)
  expect_equal(ncol(res_dense$u), k)
  expect_equal(ncol(res_dense$v), k)

  # Test sparse version
  res_sparse <- genpca:::gmd_fast_cpp(X, Qs, Rs, k)
  expect_equal(length(res_sparse$d), k)
  expect_equal(nrow(res_sparse$u), n)
  expect_equal(nrow(res_sparse$v), p)

  # Results should be very similar (same problem, different matrix types)
  expect_equal(res_dense$d, res_sparse$d, tolerance = 1e-6)

  # Check subspaces are aligned (may have sign/permutation differences)
  if (requireNamespace("clue", quietly = TRUE)) {
    align_u <- align_perm_sign(res_dense$u, res_sparse$u)
    expect_true(mean(align_u$corr) > 0.99)

    align_v <- align_perm_sign(res_dense$v, res_sparse$v)
    expect_true(mean(align_v$corr) > 0.99)
  }
})

test_that("gmd_fast_cpp handles dsyMatrix input", {
  set.seed(123)
  n <- 50
  p <- 20
  k <- 3
  X <- matrix(rnorm(n * p), n, p)

  # Create dense symmetric matrices (will become dsyMatrix)
  Q <- crossprod(matrix(rnorm(n * n), n, n)) + diag(n) * 0.1
  R <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p) * 0.1

  # Convert to Matrix format (will be dsyMatrix)
  Q_mat <- Matrix::Matrix(Q)
  R_mat <- Matrix::Matrix(R)

  expect_true(inherits(Q_mat, "dsyMatrix") || inherits(Q_mat, "dpoMatrix"))
  expect_true(inherits(R_mat, "dsyMatrix") || inherits(R_mat, "dpoMatrix"))

  # Should not error with dsyMatrix input
  res <- genpca:::gmd_fast_cpp(X, Q_mat, R_mat, k)

  expect_equal(length(res$d), k)
  expect_equal(nrow(res$u), n)
  expect_equal(nrow(res$v), p)
  expect_equal(res$k, k)
})

test_that("gmd_fast_cpp diagonal fast path matches generic path", {
  set.seed(456)
  n <- 120
  p <- 80
  k <- 6
  X <- matrix(rnorm(n * p), n, p)

  Q <- Matrix::Diagonal(n, x = runif(n, 0.8, 1.2))
  R <- Matrix::Diagonal(p, x = runif(p, 0.8, 1.2))

  res_diag <- genpca:::gmd_fast_cpp(
    X, Q, R, k,
    topk = FALSE,
    auto_topk = FALSE,
    diag_fast = TRUE
  )
  res_generic <- genpca:::gmd_fast_cpp(
    X, Q, R, k,
    topk = FALSE,
    auto_topk = FALSE,
    diag_fast = FALSE
  )

  expect_equal(res_diag$d, res_generic$d, tolerance = 1e-6)
  expect_equal(res_diag$k, res_generic$k)

  align_u <- align_perm_sign(res_diag$u, res_generic$u)
  expect_true(mean(align_u$corr) > 0.99)
  align_v <- align_perm_sign(res_diag$v, res_generic$v)
  expect_true(mean(align_v$corr) > 0.99)
})

test_that("genpca method='auto' dispatches predictably", {
  set.seed(789)

  # Small problem should choose eigen.
  X_small <- matrix(rnorm(120 * 40), 120, 40)
  fit_auto_small <- genpca::genpca(
    X_small,
    ncomp = 8,
    method = "auto",
    preproc = multivarious::center()
  )
  expect_equal(fit_auto_small$method, "eigen")

  fit_eig_small <- genpca::genpca(
    X_small,
    ncomp = 8,
    method = "eigen",
    preproc = multivarious::center()
  )
  expect_equal(fit_auto_small$sdev, fit_eig_small$sdev, tolerance = 1e-8)

  # Large diagonal, low-rank request should choose spectra when available.
  X_diag <- matrix(rnorm(260 * 220), 260, 220)
  A_diag <- Matrix::Diagonal(220, x = runif(220, 0.7, 1.3))
  M_diag <- Matrix::Diagonal(260, x = runif(260, 0.7, 1.3))
  fit_auto_diag <- genpca::genpca(
    X_diag, A = A_diag, M = M_diag,
    ncomp = 10,
    method = "auto",
    preproc = multivarious::center()
  )

  expect_equal(fit_auto_diag$method, "spectra")

  # Sparse non-diagonal constraints should prefer spectra over eigen.
  make_sparse_spd <- function(m) {
    Matrix::Diagonal(m, x = rep(1, m)) +
      Matrix::bandSparse(
        m, m,
        k = c(-1, 1),
        diagonals = list(rep(-0.2, m - 1), rep(-0.2, m - 1))
      ) +
      Matrix::Diagonal(m, x = rep(0.5, m))
  }
  X_sparse <- matrix(rnorm(150 * 900), 150, 900)
  fit_auto_sparse <- genpca::genpca(
    X_sparse,
    M = methods::as(make_sparse_spd(150), "dgCMatrix"),
    A = methods::as(make_sparse_spd(900), "dgCMatrix"),
    ncomp = 8,
    method = "auto",
    preproc = multivarious::pass()
  )
  expect_equal(fit_auto_sparse$method, "spectra")
})

test_that("genpca method='auto' selects randomized for wide sparse low-rank problems", {
  skip_if_not(exists("gmd_randomized_cpp_dn", envir = asNamespace("genpca"),
                     mode = "function"))

  set.seed(790)
  n <- 420
  p <- 4200
  k <- 12

  X <- matrix(rnorm(n * p), n, p)
  make_sparse_spd <- function(m) {
    Matrix::Diagonal(m, x = rep(1, m)) +
      Matrix::bandSparse(
        m, m,
        k = c(-1, 1),
        diagonals = list(rep(-0.2, m - 1), rep(-0.2, m - 1))
      ) +
      Matrix::Diagonal(m, x = rep(0.5, m))
  }
  M_sp <- methods::as(make_sparse_spd(n), "dgCMatrix")
  A_sp <- methods::as(make_sparse_spd(p), "dgCMatrix")

  fit_auto <- genpca::genpca(
    X,
    M = M_sp,
    A = A_sp,
    ncomp = k,
    method = "auto",
    preproc = multivarious::pass()
  )
  expect_equal(fit_auto$method, "randomized")

  fit_rand <- genpca::genpca(
    X,
    M = M_sp,
    A = A_sp,
    ncomp = k,
    method = "randomized",
    preproc = multivarious::pass()
  )
  expect_equal(fit_auto$sdev, fit_rand$sdev, tolerance = 1e-6)
})

test_that("gmd_fast_cpp dual sparse path returns metric-orthonormal factors", {
  set.seed(991)
  n <- 90
  p <- 300
  k <- 10
  X <- scale(matrix(rnorm(n * p), n, p), center = TRUE, scale = FALSE)

  make_sparse_spd <- function(m) {
    methods::as(
      Matrix::Diagonal(m, x = rep(1.5, m)) +
        Matrix::bandSparse(
          m, m,
          k = c(-1, 1),
          diagonals = list(rep(-0.2, m - 1), rep(-0.2, m - 1))
        ),
      "dgCMatrix"
    )
  }

  Q <- make_sparse_spd(n)
  R <- make_sparse_spd(p)

  res <- genpca:::gmd_fast_cpp(
    X, Q, R, k,
    topk = TRUE,
    auto_topk = TRUE,
    diag_fast = FALSE
  )

  QtU <- crossprod(res$ou, Q %*% res$ou)
  RtV <- crossprod(res$ov, R %*% res$ov)
  Ik <- diag(ncol(res$ou))

  expect_lt(max(abs(QtU - Ik)), 1e-5)
  expect_lt(max(abs(RtV - Ik)), 1e-5)
})
