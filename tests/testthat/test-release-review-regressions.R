test_that("small positive diagonal weights retain valid GMD factors", {
  X <- diag(c(1, 1e6))
  W <- diag(c(1, 1e-10))
  for (side in c("rows", "columns")) {
    M <- if (side == "rows") W else diag(2)
    A <- if (side == "columns") W else diag(2)
    for (method in c("eigen", "spectra", "randomized")) {
      # The full two-dimensional sketch needs no regularizing Gram jitter.
      f <- genpca(X, A = A, M = M, ncomp = 2, method = method, jitter_metric = 0)
      label <- paste(side, method)
      expect_equal(f$sdev, c(10, 1), tolerance = 1e-7, label = label)
      expect_equal(crossprod(f$ov, A %*% f$ov), diag(2), tolerance = 1e-7, ignore_attr = TRUE)
      expect_equal(crossprod(f$ou, M %*% f$ou), diag(2), tolerance = 1e-7, ignore_attr = TRUE)
      expect_equal(as.matrix(multivarious::reconstruct(f)), X,
                   tolerance = 1e-7, ignore_attr = TRUE)
      expect_equal(X %*% A %*% f$ov, sweep(f$ou, 2, f$sdev, "*"),
                   tolerance = 1e-7, ignore_attr = TRUE)
    }
    f <- genpca_cov(crossprod(X, M %*% X), R = A, ncomp = 2)
    expect_equal(f$d, c(10, 1), tolerance = 1e-10)
    expect_equal(as.matrix(t(f$v) %*% A %*% f$v), diag(2), tolerance = 1e-10)
    expect_equal(f$R_rank, 2L)
  }
})

test_that("diagonal PLS operators use the same support in every representation", {
  w <- c(1, 1e-10, 0)
  for (W in list(w, diag(w), Matrix::Diagonal(x = w))) {
    op <- genpca:::.metric_operators(W, n_expected = 3)
    expect_equal(as.matrix(op$mult_sqrt(op$mult_sqrt(diag(3)))), diag(w))
    expect_equal(as.matrix(op$mult_sqrt(op$mult_invsqrt(diag(3)))), diag(c(1, 1, 0)))
  }
  op <- genpca:::as_weight_operator(diag(w), sqrt = TRUE, inverse = TRUE)
  expect_equal(op(diag(3)), diag(c(1, 1e5, 0)))
  f <- gplssvd_op(diag(c(1, 1e6)), diag(2), XRW = c(1, 1e-10), k = 1)
  expect_equal(f$d, 10)
  expect_equal(as.numeric(crossprod(f$p, diag(c(1, 1e-10)) %*% f$p)), 1)
})

test_that("deflation rank decisions are independent of iteration tolerance", {
  raw <- suppressWarnings(genpca:::gmd_deflation_cpp(
    matrix(0, 5, 3), genpca:::as_dgc(diag(5)), genpca:::as_dgc(diag(3)), k = 2))
  expect_length(raw$d, 0L)
  expect_length(raw$propv, 0L)
  expect_length(raw$cumv, 0L)
  expect_equal(ncol(raw$u), 0L)
  for (cpp in c(TRUE, FALSE)) {
    set.seed(42)
    f <- suppressWarnings(genpca(diag(c(10, 1, .1)), ncomp = 3,
                                 method = "deflation", rank_rtol = .2, use_cpp = cpp))
    expect_equal(f$sdev, 10, tolerance = 1e-6)
    expect_equal(ncol(f$ov), 1L)
    expect_length(f$propv, 1L)
    expect_equal(f$cumv, cumsum(f$propv))

    for (sparse in c(FALSE, TRUE)) {
      X <- diag(c(10, 1e-7))
      if (sparse) X <- Matrix::Matrix(X, sparse = TRUE)
      set.seed(42)
      f <- genpca(X, ncomp = 2, method = "deflation", rank_rtol = 1e-10,
                   threshold = 1e-6, use_cpp = cpp)
      expect_length(f$sdev, 2L)
      expect_equal(f$sdev / c(10, 1e-7), c(1, 1), tolerance = 1e-6)
    }
  }
})

test_that("partial eigencore results cannot escape the shared backend", {
  # Inject the raw eigencore result at the RNG guard, leaving the wrapper
  # and its callers intact. No platform-dependent convergence failure needed.
  testthat::local_mocked_bindings(
    .with_rng_guard = function(expr) list(nconv = 0L), .package = "genpca")
  expect_error(genpca:::.top_svd(diag(4), 2), class = "genpca_solver_nonconvergence")
  expect_error(genpca:::.top_eigs_sym(diag(4), 2), class = "genpca_solver_nonconvergence")
  set.seed(3)
  X <- matrix(rnorm(80 * 70), 80)
  Y <- matrix(rnorm(80 * 66), 80)
  expect_error(gplssvd_op(X, Y, k = 2), class = "genpca_solver_nonconvergence")
  expect_error(genpca_cov(diag(101), ncomp = 2), class = "genpca_solver_nonconvergence")
  # Existing bounded fallbacks recover valid decompositions.
  f <- genpca:::gmd_spectra(X, Matrix::Diagonal(80), Matrix::Diagonal(70),
                            k = 2, auto_topk = FALSE)
  expect_equal(f$d, svd(X, nu = 0, nv = 0)$d[1:2], tolerance = 1e-10)
  f <- genpca:::.mnpca_init_factors(X, 2)
  expect_equal(f$X %*% t(f$W), svd(X)$u[, 1:2] %*%
                 diag(svd(X)$d[1:2]) %*% t(svd(X)$v[, 1:2]), tolerance = 1e-10)
})

test_that("explicit clipping removes negative eigenvalues within validation tolerance", {
  rot <- matrix(c(1, 1, -1, 1), 2) / sqrt(2)
  diagonal <- diag(c(1, -1e-10))
  general <- rot %*% diagonal %*% t(rot)
  for (A in list(diagonal, general, Matrix::Matrix(general, sparse = TRUE))) {
    B <- repair_metric(A, method = "clip")
    expect_true(attr(B, "repair_report")$changed)
    expect_lt(attr(B, "repair_report")$min_eigenvalue_before, 0)
    expect_gte(min(eigen(as.matrix(B), symmetric = TRUE)$values), -1e-15)
    expect_warning(f <- genpca(diag(2), A = A, ncomp = 1, constraints_remedy = "clip"),
                   class = "genpca_metric_repaired")
    expect_gte(min(eigen(as.matrix(f$A), symmetric = TRUE)$values), -1e-15)
    expect_warning(genpca:::.metric_operators(A, remedy = "clip"),
                   class = "genpca_metric_repaired")
  }
  B <- repair_metric(diag(c(1, 0)), method = "clip")
  expect_false(attr(B, "repair_report")$changed)
})

test_that("geigen_cov solves the projected equation for singular metrics", {
  C0 <- matrix(c(2, 1, 1, 2), 2)
  R0 <- diag(c(1, 0))
  for (U in list(diag(2), matrix(c(1, 1, -1, 1), 2) / sqrt(2))) {
    C <- U %*% C0 %*% t(U)
    R <- U %*% R0 %*% t(U)
    f <- geigen_cov(C, R, ncomp = 1)
    residual <- C %*% f$v - (R %*% f$v) * f$lambda
    expect_equal(as.numeric(R %*% residual), c(0, 0), tolerance = 1e-12)
    expect_equal(as.numeric(R %*% f$v), as.numeric(f$v), tolerance = 1e-12)
    expect_equal(as.numeric(t(f$v) %*% R %*% f$v), 1, tolerance = 1e-12)
    expect_equal(sum(residual^2), 1, tolerance = 1e-12)
  }
  C <- diag(c(4, 1)); R <- diag(c(9, 1))
  gmd <- genpca_cov(C, R, ncomp = 1)
  ge <- geigen_cov(C, R, ncomp = 1)
  expect_equal(which.max(abs(as.numeric(gmd$v))), 1L)
  expect_equal(which.max(abs(as.numeric(ge$v))), 2L)
  expect_equal(as.numeric(C %*% ge$v - (R %*% ge$v) * ge$lambda), c(0, 0))
})
