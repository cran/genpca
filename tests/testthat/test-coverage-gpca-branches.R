test_that("prep_constraints covers matrix validation, identity remedy, and preprocessor detection", {
  X <- matrix(rnorm(20), 5, 4)

  nonsym_A <- diag(4)
  nonsym_A[1, 2] <- 0.25
  expect_error(genpca:::prep_constraints(X, nonsym_A, NULL), "Matrix A must be symmetric")

  nonsym_M <- diag(5)
  nonsym_M[1, 2] <- 0.25
  expect_error(genpca:::prep_constraints(X, NULL, nonsym_M), "Matrix M must be symmetric")

  bad_A <- diag(c(1, -1, 1, 1))
  bad_M <- diag(c(1, -1, 1, 1, 1))
  expect_message(
    pc <- genpca:::prep_constraints(
      X,
      A = bad_A,
      M = bad_M,
      remedy = "identity",
      verbose = TRUE
    ),
    "Matrix A is not SPD"
  )
  expect_s4_class(pc$A, "ddiMatrix")
  expect_s4_class(pc$M, "ddiMatrix")

  sym_A <- Matrix::Matrix(crossprod(matrix(rnorm(16), 4, 4)) + diag(4), sparse = FALSE)
  sym_M <- Matrix::Matrix(crossprod(matrix(rnorm(25), 5, 5)) + diag(5), sparse = FALSE)
  pc_dense <- genpca:::prep_constraints(X, sym_A, sym_M)
  expect_s4_class(pc_dense$A, "dgeMatrix")
  expect_s4_class(pc_dense$M, "dgeMatrix")

  pp <- multivarious::pass()
  expect_true(genpca:::is_pass_preproc(pp))
  ft <- multivarious::fit_transform(pp, X)
  expect_true(genpca:::is_pass_preproc(ft$preproc))
})

test_that("genpca covers sparse pass-through, verbose, deflation, spectra, randomized, truncation, and reconstruction branches", {
  set.seed(901)
  X <- matrix(rnorm(8 * 5), 8, 5)
  rownames(X) <- paste0("obs", seq_len(nrow(X)))
  colnames(X) <- paste0("var", seq_len(ncol(X)))

  expect_message(
    fit_verbose <- genpca::genpca(
      X,
      ncomp = 2,
      method = "eigen",
      preproc = multivarious::pass(),
      verbose = TRUE
    ),
    "GPCA finished"
  )
  expect_equal(fit_verbose$method, "eigen")
  expect_equal(rownames(fit_verbose$v), colnames(X))
  expect_error(genpca:::truncate.genpca(fit_verbose), "ncomp")
  expect_error(genpca:::truncate.genpca(fit_verbose, ncomp = 99), "positive integer")
  expect_identical(genpca:::truncate.genpca(fit_verbose, ncomp = multivarious::ncomp(fit_verbose)), fit_verbose)
  fit_one <- genpca:::truncate.genpca(fit_verbose, ncomp = 1)
  expect_equal(multivarious::ncomp(fit_one), 1L)
  expect_error(multivarious::reconstruct(fit_verbose, comp = 99), "out of bounds")
  rec_subset <- suppressWarnings(multivarious::reconstruct(fit_verbose, comp = 1, rowind = 1:3, colind = 1:2))
  expect_equal(dim(rec_subset), c(3L, 2L))

  X_sp <- methods::as(Matrix::Matrix(X, sparse = TRUE), "dgCMatrix")
  fit_sparse <- genpca::genpca(
    X_sp,
    ncomp = 2,
    method = "eigen",
    preproc = multivarious::pass()
  )
  expect_equal(fit_sparse$method, "eigen")

  fit_defl_primal <- genpca::genpca(
    X,
    ncomp = 2,
    method = "deflation",
    use_cpp = TRUE,
    maxit_deflation = 50,
    preproc = multivarious::pass()
  )
  expect_equal(fit_defl_primal$method, "deflation")

  X_wide <- matrix(rnorm(5 * 8), 5, 8)
  fit_defl_dual <- genpca::genpca(
    X_wide,
    ncomp = 2,
    method = "deflation",
    use_cpp = TRUE,
    maxit_deflation = 50,
    preproc = multivarious::pass()
  )
  expect_equal(fit_defl_dual$method, "deflation")

  fit_rand <- genpca::genpca(
    X,
    ncomp = 2,
    method = "randomized",
    oversample = 2,
    n_power = 1,
    n_polish = 1,
    seed_randomized = 11,
    preproc = multivarious::pass(),
    verbose = TRUE
  )
  expect_equal(fit_rand$method, "randomized")

  fit_zero <- suppressWarnings(genpca::genpca(
    matrix(0, 5, 4),
    ncomp = 2,
    method = "spectra",
    preproc = multivarious::pass()
  ))
  expect_equal(multivarious::ncomp(fit_zero), 0L)
})

test_that("gmdLA covers diagonal errors, non-diagonal primal and dual paths, approximations, and zero variance warnings", {
  set.seed(902)
  X <- matrix(rnorm(8 * 5), 8, 5)
  Q <- Matrix::Diagonal(8)
  R_bad <- Matrix::Diagonal(5, x = c(1, 1, 1, 1, -1))
  expect_error(
    genpca:::gmdLA(X, Q, R_bad, k = 2, n_orig = 8, p_orig = 5, use_dual = FALSE),
    "must be PSD|non-negative"
  )

  make_spd <- function(m) {
    B <- matrix(rnorm(m * m), m, m)
    Matrix::Matrix(crossprod(B) + diag(m), sparse = FALSE)
  }
  R <- make_spd(5)
  fit_primal <- genpca:::gmdLA(
    X, Q, R,
    k = 2,
    n_orig = 8,
    p_orig = 5,
    use_dual = FALSE,
    verbose = TRUE
  )
  expect_equal(fit_primal$k, 2L)

  X_dual <- matrix(rnorm(5 * 8), 5, 8)
  Q_dual <- make_spd(5)
  R_dual <- Matrix::Diagonal(8)
  fit_dual <- genpca:::gmdLA(
    X_dual, Q_dual, R_dual,
    k = 2,
    n_orig = 5,
    p_orig = 8,
    use_dual = TRUE,
    verbose = TRUE
  )
  expect_equal(fit_dual$k, 2L)

  # A positive definite metric larger than maxeig is factored exactly
  # (Cholesky), never truncated: the result matches the unrestricted call.
  X_big <- matrix(rnorm(130 * 105), 130, 105)
  R_big_raw <- matrix(rnorm(105 * 6), 105, 6)
  R_big <- Matrix::Matrix(crossprod(t(R_big_raw)) + diag(105), sparse = FALSE)
  fit_exact <- genpca:::gmdLA(
    X_big, Matrix::Diagonal(130), R_big,
    k = 2, n_orig = 130, p_orig = 105, maxeig = Inf, use_dual = FALSE
  )
  expect_silent(
    fit_guard <- genpca:::gmdLA(
      X_big, Matrix::Diagonal(130), R_big,
      k = 2, n_orig = 130, p_orig = 105, maxeig = 20, use_dual = FALSE
    )
  )
  expect_equal(fit_guard$k, 2L)
  expect_equal(fit_guard$d, fit_exact$d, tolerance = 1e-8)
  # A singular general metric above maxeig needs a dense eigendecomposition
  # and is refused rather than approximated.
  R_sing <- Matrix::Matrix(crossprod(t(R_big_raw)), sparse = FALSE)
  expect_error(
    genpca:::gmdLA(
      X_big, Matrix::Diagonal(130), R_sing,
      k = 2, n_orig = 130, p_orig = 105, maxeig = 20, use_dual = FALSE
    ),
    "maxeig"
  )
  fit_sing <- genpca:::gmdLA(
    X_big, Matrix::Diagonal(130), R_sing,
    k = 2, n_orig = 130, p_orig = 105, maxeig = Inf, use_dual = FALSE
  )
  expect_equal(fit_sing$k, 2L)

  expect_error(
    genpca:::gmdLA(
      matrix(0, 4, 3),
      Matrix::Diagonal(4),
      Matrix::Diagonal(3),
      k = 2,
      n_orig = 4,
      p_orig = 3,
      use_dual = FALSE
    ),
    "No positive eigenvalues"
  )
})

test_that("gpca_mle covers determinant scaling and convergence messages", {
  set.seed(903)
  X <- matrix(rnorm(6 * 4), 6, 4)
  fit <- suppressWarnings(suppressMessages(genpca::gpca_mle(
      X,
      ncomp = 2,
      max_iter = 2,
      scale_fix = "det",
      tol = Inf,
      method = "eigen",
      preproc = multivarious::pass(),
      verbose = TRUE
    )))
  expect_length(fit$loglik, 1L)
  expect_lte(length(fit$loglik_path), 2L)
})

test_that("genpca_cov covers GMD and generalized-eigen covariance branches", {
  set.seed(904)
  X <- matrix(rnorm(20 * 6), 20, 6)
  C <- crossprod(scale(X, center = TRUE, scale = FALSE))
  R_general <- crossprod(matrix(rnorm(36), 6, 6)) + diag(6)
  R_nonsym <- R_general
  R_nonsym[1, 2] <- R_nonsym[1, 2] + 0.5

  expect_error(genpca::genpca_cov(C[1:5, ], ncomp = 2), "p == ncol")
  expect_error(genpca::genpca_cov(C, R = c(1, 1), ncomp = 2), "length")
  expect_error(genpca::genpca_cov(C, R = c(1, -2, rep(1, 4)), ncomp = 2), "non-negative")

  # Asymmetric metrics are an input error, not silently symmetrized.
  expect_error(genpca::genpca_cov(C, R = R_nonsym, ncomp = 3), "symmetric")
  fit_gmd <- genpca::genpca_cov(C, R = R_general, ncomp = 3, verbose = TRUE)
  expect_equal(fit_gmd$method, "gmd")
  expect_equal(fit_gmd$k, 3L)
  expect_equal(as.matrix(crossprod(fit_gmd$v, R_general %*% fit_gmd$v)),
               diag(3), tolerance = 1e-6)

  C_big <- diag(seq(120, 1, length.out = 120))
  fit_big <- genpca::genpca_cov(C_big, R = NULL, ncomp = 2)
  expect_equal(fit_big$k, 2L)

  expect_warning(
    fit_warn <- genpca::geigen_cov(
      matrix(c(1, 0, 0, -1), 2, 2),
      R = NULL,
      ncomp = 1
    ),
    "non-PSD"
  )
  expect_warning(
    genpca::genpca_cov(C, R = NULL, ncomp = 1, method = "geigen"),
    "deprecated"
  )
  expect_equal(fit_warn$method, "geigen")

  R_bad_diag <- Matrix::Diagonal(6, x = c(1, -0.5, rep(1, 4)))
  expect_error(
    genpca::geigen_cov(C, R = R_bad_diag, ncomp = 2, constraints_remedy = "error"),
    "positive semi-definite|PSD"
  )
  fit_ridge <- genpca::geigen_cov(C, R = R_bad_diag, ncomp = 2, constraints_remedy = "ridge", verbose = TRUE)
  fit_clip <- genpca::geigen_cov(C, R = R_bad_diag, ncomp = 2, constraints_remedy = "clip", verbose = TRUE)
  fit_identity <- genpca::geigen_cov(C, R = R_bad_diag, ncomp = 2, constraints_remedy = "identity", verbose = TRUE)
  expect_equal(fit_ridge$k, 2L)
  expect_equal(fit_clip$k, 2L)
  expect_equal(fit_identity$k, 2L)

  R_indef <- matrix(c(1, 2, 2, 1), 2, 2)
  C2 <- diag(2)
  fit_general_clip <- genpca::geigen_cov(
    C2,
    R = R_indef,
    ncomp = 1,
    constraints_remedy = "clip",
    verbose = TRUE
  )
  fit_general_identity <- genpca::geigen_cov(
    C2,
    R = R_indef,
    ncomp = 1,
    constraints_remedy = "identity",
    verbose = TRUE
  )
  expect_equal(fit_general_clip$k, 1L)
  expect_equal(fit_general_identity$R_rank, 2L)
})

test_that("C++ GMD top-k and randomized branches cover operator and zero-output paths", {
  set.seed(905)
  make_spd <- function(m) {
    B <- matrix(rnorm(m * m), m, m)
    crossprod(B) + diag(m)
  }
  make_sparse_spd <- function(m) {
    methods::as(
      Matrix::Diagonal(m, x = rep(1.5, m)) +
        Matrix::bandSparse(m, m, k = c(-1, 1),
                           diagonals = list(rep(-0.2, m - 1), rep(-0.2, m - 1))),
      "dgCMatrix"
    )
  }

  X_primal <- matrix(rnorm(80 * 30), 80, 30)
  res_top_primal <- genpca:::gmd_fast_cpp(
    X_primal,
    Matrix::Matrix(make_spd(80), sparse = FALSE),
    Matrix::Matrix(make_spd(30), sparse = FALSE),
    k = 3,
    topk = TRUE,
    auto_topk = FALSE,
    diag_fast = FALSE,
    cache = FALSE
  )
  expect_equal(res_top_primal$k, 3L)

  X_dual <- matrix(rnorm(30 * 80), 30, 80)
  res_top_dual <- genpca:::gmd_fast_cpp(
    X_dual,
    Matrix::Matrix(make_spd(30), sparse = FALSE),
    Matrix::Matrix(make_spd(80), sparse = FALSE),
    k = 3,
    topk = TRUE,
    auto_topk = FALSE,
    diag_fast = FALSE,
    cache = FALSE
  )
  expect_equal(res_top_dual$k, 3L)

  zero_rand <- genpca:::gmd_randomized_cpp_dn(
    matrix(0, 7, 5),
    diag(7),
    diag(5),
    k = 2,
    oversample = 1,
    n_power = 1
  )
  expect_equal(zero_rand$k, 0L)

  X <- matrix(rnorm(7 * 5), 7, 5)
  Qs <- make_sparse_spd(7)
  Rs <- make_sparse_spd(5)
  Qd <- Matrix::Matrix(make_spd(7), sparse = FALSE)
  Rd <- Matrix::Matrix(make_spd(5), sparse = FALSE)

  res_q_sparse <- genpca:::gmd_randomized(X, Qs, Rd, k = 2, oversample = 1, n_power = 1, seed = 13, use_cpp = TRUE)
  res_r_sparse <- genpca:::gmd_randomized(X, Qd, Rs, k = 2, oversample = 1, n_power = 1, seed = 13, use_cpp = TRUE)
  res_both_sparse <- genpca:::gmd_randomized(X, Qs, Rs, k = 2, oversample = 1, n_power = 1, seed = 13, use_cpp = TRUE)
  expect_equal(res_q_sparse$k, 2L)
  expect_equal(res_r_sparse$k, 2L)
  expect_equal(res_both_sparse$k, 2L)
})
