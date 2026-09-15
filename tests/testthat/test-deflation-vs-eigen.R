testthat::test_that("deflation closely matches eigen on small dense problems", {
  skip_if_not_installed("Matrix")
  library(Matrix)

  set.seed(123)
  n <- 60
  p <- 40
  k <- 8
  X <- matrix(rnorm(n * p), n, p)

  # Eigen (reference) vs Deflation (R path)
  fit_eig <- genpca(X, ncomp = k, method = "eigen",
                    preproc = multivarious::center())
  fit_def <- genpca(X, ncomp = k, method = "deflation", use_cpp = FALSE,
                    preproc = multivarious::center(), threshold = 1e-8)

  # 1) Singular values: tight numerical match
  testthat::expect_equal(fit_eig$sdev[1:k], fit_def$sdev[1:k], tolerance = 1e-4)

  # 2) Subspace agreement via principal angles (scores space)
  # Use ou (orthonormal in M metric). With identity metrics here, it is Euclidean-orthonormal.
  Ou_e <- as.matrix(fit_eig$ou)[, 1:k, drop = FALSE]
  Ou_d <- as.matrix(fit_def$ou)[, 1:k, drop = FALSE]
  sv  <- svd(t(Ou_e) %*% Ou_d)$d
  # Sum of squared canonical correlations should be ~ k
  testthat::expect_equal(sum(sv^2), k, tolerance = 1e-3)

  # 3) Scores and loadings: allow arbitrary orthogonal rotation within the subspace.
  procrustes_align <- function(B, A) {
    # Find orthogonal R minimizing ||A - B R||_F via SVD(B^T A).
    BtA <- crossprod(B, A)
    sv  <- svd(BtA)
    R   <- sv$u %*% t(sv$v)
    B %*% R
  }

  S_e <- as.matrix(multivarious::scores(fit_eig))[, 1:k, drop = FALSE]
  S_d <- as.matrix(multivarious::scores(fit_def))[, 1:k, drop = FALSE]
  S_d_al <- procrustes_align(S_d, S_e)
  testthat::expect_lt(norm(S_e - S_d_al, "F"), 1e-3 * (1 + norm(S_e, "F")))

  L_e <- as.matrix(multivarious::components(fit_eig))[, 1:k, drop = FALSE]
  L_d <- as.matrix(multivarious::components(fit_def))[, 1:k, drop = FALSE]
  L_d_al <- procrustes_align(L_d, L_e)
  testthat::expect_lt(norm(L_e - L_d_al, "F"), 1e-3 * (1 + norm(L_e, "F")))
})

testthat::test_that("C++ deflation closely matches eigen on small dense problems", {
  skip_if_not_installed("Matrix")

  set.seed(321)
  n <- 55
  p <- 35
  k <- 6
  X <- matrix(rnorm(n * p), n, p)

  fit_eig <- genpca(X, ncomp = k, method = "eigen",
                    preproc = multivarious::center())
  fit_cpp <- genpca(X, ncomp = k, method = "deflation", use_cpp = TRUE,
                    preproc = multivarious::center(), threshold = 1e-8,
                    maxit_deflation = 700L)

  testthat::expect_equal(fit_eig$sdev[1:k], fit_cpp$sdev[1:k], tolerance = 1e-4)
})

testthat::test_that("C++ deflation dispatch keeps sparse X sparse", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("testthat", minimum_version = "3.0.0")
  library(Matrix)

  X <- rsparsematrix(12, 9, density = 0.2)
  Q <- genpca:::as_dgc(Diagonal(nrow(X)))
  R <- genpca:::as_dgc(Diagonal(ncol(X)))
  saw_sparse <- FALSE

  testthat::local_mocked_bindings(
    gmd_deflation_cpp_sp = function(X, Q, R, k, thr = 1e-7,
                                    maxit = 500L, verbose = FALSE, rank_rtol = 1e-6) {
      saw_sparse <<- methods::is(X, "sparseMatrix")
      list(d = numeric(0), v = matrix(0, ncol(X), 0),
           u = matrix(0, nrow(X), 0), k = 0,
           cumv = numeric(0), propv = numeric(0))
    },
    .package = "genpca"
  )

  res <- genpca:::gmd_deflation_cpp_dispatch(X, Q, R, k = 2,
                                             thr = 1e-8, maxit = 20L,
                                             verbose = FALSE)

  testthat::expect_true(saw_sparse)
  testthat::expect_equal(res$k, 0)
})

testthat::test_that("sparse C++ deflation path agrees with dense C++ path", {
  skip_if_not_installed("Matrix")
  library(Matrix)

  set.seed(456)
  X_sparse <- rsparsematrix(45, 30, density = 0.12)
  X_dense <- as.matrix(X_sparse)
  Q <- Diagonal(nrow(X_sparse))
  R <- Diagonal(ncol(X_sparse))
  k <- 4

  res_sparse <- genpca:::gmd_deflation_cpp_sp(
    genpca:::as_dgc(X_sparse), genpca:::as_dgc(Q), genpca:::as_dgc(R),
    k = k, thr = 1e-8, maxit = 700L
  )
  res_dense <- genpca:::gmd_deflation_cpp(
    X_dense, genpca:::as_dgc(Q), genpca:::as_dgc(R),
    k = k, thr = 1e-8, maxit = 700L
  )

  testthat::expect_equal(res_sparse$d, res_dense$d, tolerance = 1e-5)
})

testthat::test_that("genpca C++ deflation accepts sparse pass-through input", {
  skip_if_not_installed("Matrix")
  library(Matrix)

  set.seed(789)
  X_sparse <- rsparsematrix(40, 25, density = 0.1)
  k <- 3

  fit_sparse <- genpca(
    X_sparse, ncomp = k, method = "deflation", use_cpp = TRUE,
    preproc = multivarious::pass(), threshold = 1e-8,
    maxit_deflation = 700L
  )
  fit_dense <- genpca(
    as.matrix(X_sparse), ncomp = k, method = "deflation", use_cpp = TRUE,
    preproc = multivarious::pass(), threshold = 1e-8,
    maxit_deflation = 700L
  )

  testthat::expect_equal(fit_sparse$sdev, fit_dense$sdev, tolerance = 1e-5)
  testthat::expect_equal(dim(multivarious::scores(fit_sparse)), c(nrow(X_sparse), k))
})
