# Behaviours the eigencore wrapper must provide (plans/2026-09-04-eigencore-swap.md, section 3)

test_that(".top_svd and .top_eigs_sym never disturb .Random.seed", {
  X <- matrix(rnorm(400 * 30), 400)
  set.seed(11); before <- .Random.seed
  invisible(genpca:::.top_svd(X, 3))
  invisible(genpca:::.top_eigs_sym(crossprod(X), 3, "LA"))
  invisible(genpca:::.top_svd(function(x, args) X %*% x, 2,
                              adjoint = function(x, args) crossprod(X, x), dim = dim(X)))
  expect_identical(.Random.seed, before)
  set.seed(11); r1 <- runif(1)
  set.seed(11); invisible(genpca:::.top_svd(X, 3)); r2 <- runif(1)
  expect_identical(r1, r2)
})

test_that(".top_svd wraps vector-returning operators and passes args", {
  set.seed(2); X <- matrix(rnorm(200 * 20), 200)
  ref <- svd(X)
  sv <- genpca:::.top_svd(function(x, args) as.numeric(args$scale * (X %*% x)), 3,
                          adjoint = function(x, args) as.numeric(args$scale * crossprod(X, x)),
                          dim = dim(X), args = list(scale = 2))
  expect_equal(sv$d, 2 * ref$d[1:3], tolerance = 1e-8)
  expect_true(sv$converged)
  expect_equal(dim(sv$u), c(200L, 3L))
  expect_equal(dim(sv$v), c(20L, 3L))
})

test_that(".top_eigs_sym honours tol and which, and reports convergence", {
  set.seed(3); S <- crossprod(matrix(rnorm(300 * 150), 300)) - 0.5 * diag(150)
  ref <- eigen(S, symmetric = TRUE)
  la <- genpca:::.top_eigs_sym(S, 4, "LA", tol = 1e-12)
  expect_equal(la$values, ref$values[1:4], tolerance = 1e-10)
  expect_true(la$converged)
  sa <- genpca:::.top_eigs_sym(S, 1, "SA")
  expect_equal(sa$values, min(ref$values), tolerance = 1e-8)
  lm <- genpca:::.top_eigs_sym(S, 2, "LM")
  expect_equal(sort(abs(lm$values), decreasing = TRUE), sort(abs(ref$values), decreasing = TRUE)[1:2], tolerance = 1e-8)
})

test_that("gplssvd_op forwards svd_opts$tol to the eigencore backend", {
  set.seed(4); n <- 150; p <- 80; q <- 70
  X <- matrix(rnorm(n * p), n); Y <- matrix(rnorm(n * q), n)
  loose <- gplssvd_op(X, Y, k = 3, svd_opts = list(tol = 1e-3))
  tight <- gplssvd_op(X, Y, k = 3, svd_opts = list(tol = 1e-12))
  ref <- svd(crossprod(X, Y))$d[1:3]
  expect_equal(tight$d, ref, tolerance = 1e-10)
  expect_equal(loose$d, ref, tolerance = 1e-2)
  expect_equal(gplssvd_op(X, Y, k = 3, svd_backend = "RSpectra")$d, tight$d, tolerance = 1e-8)
})

test_that("the whitened GMD operator converges under eigencore and matches the dense SVD", {
  set.seed(5); n <- 400; p <- 60; k <- 5
  X <- matrix(rnorm(n * p), n)
  R <- crossprod(matrix(rnorm(p * p), p)) / p + 0.1 * diag(p)
  Qd <- runif(n, 0.5, 2)
  FR <- genpca:::.metric_factor(R); FQ <- genpca:::.metric_factor(Matrix::Diagonal(n, x = Qd))
  op <- function(V, args = NULL) FQ$apply_t(X %*% FR$apply(V))
  opt <- function(U, args = NULL) FR$apply_t(crossprod(X, FQ$apply(U)))
  sv <- genpca:::.top_svd(op, k, tol = 1e-10, adjoint = opt, dim = c(FQ$ncol, FR$ncol))
  expect_true(sv$converged)
  expect_equal(sv$nconv, k)
  B <- as.matrix(FQ$apply_t(X %*% FR$mat))
  expect_equal(sv$d, svd(B)$d[1:k], tolerance = 1e-9)
})

test_that("spectra never factors a singular large-side metric and matches eigen", {
  set.seed(6); n <- 300; p <- 20; k <- 3
  X <- matrix(rnorm(n * p), n)
  # graph Laplacian on the rows: PSD, exactly singular, on the large side
  W <- Matrix::bandSparse(n, n, k = 1, diagonals = list(rep(1, n - 1)), symmetric = TRUE)
  L <- Matrix::Diagonal(n, x = Matrix::rowSums(W)) - W
  f_sp <- genpca(X, M = L, ncomp = k, method = "spectra", preproc = multivarious::pass(), maxeig = 50)
  f_ei <- genpca(X, M = L, ncomp = k, method = "eigen", preproc = multivarious::pass())
  expect_equal(f_sp$sdev, f_ei$sdev, tolerance = 1e-6)
  # a singular metric on the SMALL side above maxeig would need a dense
  # eigendecomposition: refused with guidance under eigen and spectra, and
  # routed to deflation by auto
  Xs <- matrix(rnorm(300 * 100), 300)
  W2 <- Matrix::bandSparse(100, 100, k = 1, diagonals = list(rep(1, 99)), symmetric = TRUE)
  L2 <- Matrix::Diagonal(100, x = Matrix::rowSums(W2)) - W2
  expect_error(genpca(Xs, A = L2, ncomp = k, method = "spectra", preproc = multivarious::pass(), maxeig = 50), "maxeig")
  expect_error(genpca(Xs, A = L2, ncomp = k, method = "eigen", preproc = multivarious::pass(), maxeig = 50), "maxeig")
  f_auto <- suppressWarnings(genpca(Xs, A = L2, ncomp = k, method = "auto", preproc = multivarious::pass(), maxeig = 50))
  expect_equal(f_auto$method, "deflation")
  f_ok <- genpca(Xs, A = L2, ncomp = k, method = "spectra", preproc = multivarious::pass(), maxeig = Inf)
  expect_equal(f_ok$sdev, genpca(Xs, A = L2, ncomp = k, method = "eigen", preproc = multivarious::pass())$sdev, tolerance = 1e-6)
})
