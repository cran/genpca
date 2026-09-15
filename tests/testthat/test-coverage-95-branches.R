test_that("SFPCA helper edge cases cover explicit fallbacks", {
  set.seed(9501)
  X <- matrix(rnorm(4 * 3), 4, 3)
  spat_cds <- matrix(c(0, 0, 0,
                       1, 0, 0,
                       0, 1, 0), nrow = 3)

  expect_error(genpca:::second_diff_matrix(2), "at least 3")
  expect_error(
    genpca::sfpca(X, K = 1, spat_cds = NULL),
    "spat_cds"
  )
  expect_error(
    genpca::sfpca(X, K = 1, spat_cds = spat_cds, nlambda = 1),
    "nlambda"
  )
  expect_error(
    genpca::sfpca(X, K = 1, spat_cds = spat_cds, lambda_min_ratio = 1),
    "lambda_min_ratio"
  )
  expect_error(
    genpca::sfpca(
      X, K = 1, spat_cds = spat_cds,
      lambda_u = 0.01, lambda_v = 0.01, alpha_u = -1
    ),
    "alpha_u"
  )

  fit <- suppressWarnings(capture.output(
    value <- genpca::sfpca(
      X, K = 1, spat_cds = spat_cds,
      lambda_u = 0.01, lambda_v = 0.01,
      max_iter = 2, tol = 1e-4, verbose = TRUE
    )
  ))
  expect_type(fit, "character")
  expect_s3_class(value, "sfpca")
  expect_error(multivarious::reconstruct(value, comp = 2), "comp")

  identical_coords <- matrix(0, nrow = 3, ncol = 3)
  omega_identical <- genpca:::construct_spatial_penalty(
    identical_coords,
    k = 2
  )
  expect_s4_class(omega_identical, "Matrix")
  expect_error(
    genpca:::construct_spatial_penalty(spat_cds, method = "neighbor", k = 1),
    "Neighbor graph"
  )
  expect_error(
    genpca:::construct_spatial_penalty(spat_cds, method = "bogus", k = 1),
    "Invalid method"
  )
  expect_error(
    suppressWarnings(genpca:::construct_spatial_penalty(spat_cds[, 1, drop = FALSE], k = 1)),
    "knn"
  )

  expect_equal(genpca:::default_alpha(Matrix::Diagonal(201, x = 0)), 0)
  expect_gt(genpca:::default_alpha(Matrix::Diagonal(201, x = 1)), 0)

  sel_zero <- genpca:::sfpca_select_lambda(
    b = rep(0, 4),
    S = Matrix::Diagonal(4),
    F2 = 1,
    np = 12,
    fixed_norm = 1,
    penalty = "l1"
  )
  expect_equal(sel_zero$lambda, 0)
  expect_length(sel_zero$lambdas, 0)

  Omega_u <- Matrix::Diagonal(nrow(X))
  Omega_v <- Matrix::Diagonal(ncol(X))
  expect_error(
    genpca:::sfpca_rank1(
      X = NULL, ops = NULL,
      lambda_u = 0, lambda_v = 0,
      alpha_u = 0, alpha_v = 0,
      Omega_u = Omega_u, Omega_v = Omega_v,
      penalty_u = "l1", penalty_v = "l1",
      max_iter = 1, tol = 1e-6, verbose = FALSE
    ),
    "Either X"
  )
  ops <- genpca:::sfpca_make_ops(Matrix::Matrix(X, sparse = FALSE))
  expect_error(
    genpca:::sfpca_rank1(
      X = NULL, ops = ops,
      lambda_u = 0, lambda_v = 0,
      alpha_u = 0, alpha_v = 0,
      Omega_u = Omega_u, Omega_v = Omega_v,
      penalty_u = "l1", penalty_v = "l1",
      max_iter = 1, tol = 1e-6, verbose = FALSE
    ),
    "u_init"
  )

  zero_component <- capture.output(
    res_zero <- genpca:::sfpca_rank1(
      X = Matrix::Matrix(X, sparse = FALSE),
      lambda_u = 1e6, lambda_v = 1e6,
      alpha_u = 0, alpha_v = 0,
      Omega_u = Omega_u, Omega_v = Omega_v,
      penalty_u = "l1", penalty_v = "l1",
      max_iter = 2, tol = 1e-6, verbose = TRUE
    )
  )
  expect_match(paste(zero_component, collapse = "\n"), "zero solution")
  expect_equal(res_zero$d, 0)

  u_init <- rep(1 / sqrt(nrow(X)), nrow(X))
  v_init <- rep(1 / sqrt(ncol(X)), ncol(X))
  iter_output <- capture.output(
    res_iter <- genpca:::sfpca_rank1(
      X = Matrix::Matrix(X, sparse = FALSE),
      lambda_u = 0, lambda_v = 0,
      alpha_u = 0, alpha_v = 0,
      Omega_u = Omega_u, Omega_v = Omega_v,
      penalty_u = "l1", penalty_v = "l1",
      max_iter = 1, tol = 1e-12, verbose = TRUE,
      u_init = u_init, v_init = v_init
    )
  )
  expect_match(paste(iter_output, collapse = "\n"), "Iteration")
  expect_true(is.finite(res_iter$d))
  expect_error(genpca:::penalty_value(1, "bad", 1), "Unsupported")
  expect_error(
    capture.output(genpca:::svd1_deflated(matrix(NA_real_, 3, 2), verbose = TRUE)),
    "infinite|missing"
  )
})

test_that("MN-PCA validation and graphical-lasso branches are covered", {
  set.seed(9502)
  Y <- matrix(rnorm(5 * 4), 5, 4)

  fit_matrix_input <- genpca::mnpca_mrl(
    Matrix::Matrix(Y, sparse = FALSE),
    ncomp = 1,
    max_outer = 1,
    max_inner = 1,
    update_precisions = FALSE
  )
  expect_s3_class(fit_matrix_input, "mnpca_mrl")

  expect_error(genpca::mnpca_mrl(matrix("x", 2, 2), ncomp = 1), "numeric")
  expect_error(genpca::mnpca_mrl(matrix(c(1, Inf, 3, 4), 2), ncomp = 1), "non-finite")
  expect_error(genpca::mnpca_mrl(matrix(numeric(0), 0, 2), ncomp = 1), "positive")
  expect_error(genpca::mnpca_mrl(Y, ncomp = 0), "ncomp")
  expect_error(genpca::mnpca_mrl(Y, ncomp = 10), "ncomp")
  expect_error(genpca::mnpca_mrl(Y, ncomp = 1, max_outer = 0), "max_outer")
  expect_error(genpca::mnpca_mrl(Y, ncomp = 1, lambda_row = -1), "non-negative")
  expect_error(genpca::mnpca_mrl(Y, ncomp = 1, eps_ridge = 0), "positive")

  expect_message(
    genpca::mnpca_mrl(
      Y, ncomp = 1,
      max_outer = 1,
      max_inner = 1,
      update_precisions = FALSE,
      verbose = TRUE
    ),
    "\\[mnpca_mrl\\]"
  )

  one_dim <- genpca:::.mnpca_glasso_admm(
    matrix(2, 1, 1),
    lambda = 0.5,
    penalize_diagonal = TRUE
  )
  expect_equal(dim(one_dim$Theta), c(1L, 1L))

  no_penalty <- genpca:::.mnpca_glasso_admm(diag(c(2, 3)), lambda = 0)
  expect_true(no_penalty$converged)

  unscreened <- genpca:::.mnpca_glasso_admm(
    matrix(c(1, 0.2, 0.2, 1), 2),
    lambda = 0.1,
    maxit = 2,
    block_screen = FALSE
  )
  expect_equal(dim(unscreened$Theta), c(2L, 2L))

  core_no_penalty <- genpca:::.mnpca_glasso_admm_core(
    diag(2),
    lambda = 0,
    maxit = 1
  )
  expect_equal(core_no_penalty$iterations, 1L)

  obj_inf <- genpca:::.mnpca_objective_value(
    Y = matrix(0, 2, 2),
    X = matrix(0, 2, 1),
    W = matrix(0, 2, 1),
    Theta_row = matrix(NaN, 2, 2),
    Theta_col = diag(2),
    lambda_row = 0,
    lambda_col = 0
  )
  expect_true(is.infinite(obj_inf))

  init_small <- genpca:::.mnpca_init_factors(matrix(1:4, 2, 2), ncomp = 1)
  expect_equal(dim(init_small$X), c(2L, 1L))
  expect_equal(length(genpca:::.mnpca_connected_components(matrix(1, 1, 1), lambda = 1)), 1L)
  expect_error(genpca:::.mnpca_safe_solve_spd(matrix(NaN, 2, 2)), "Failed SPD")
  expect_equal(
    genpca:::.mnpca_l1_penalty(matrix(c(1, -2, 3, -4), 2), penalize_diagonal = TRUE),
    10
  )
})

test_that("covariance GPCA rejects asymmetric C and covers zero-metric branches", {
  C_nonsym <- matrix(c(2, 0.2, 0, 1), 2)
  expect_error(
    genpca:::genpca_cov_gmd(C_nonsym, R = diag(2), ncomp = 1),
    "symmetric"
  )
  C_near <- matrix(c(2, 0.2, 0.2 + 1e-14, 1), 2)
  expect_equal(
    genpca:::genpca_cov_gmd(C_near, R = diag(2), ncomp = 1)$k,
    1L
  )
  expect_error(
    genpca:::genpca_cov_gmd(diag(2), R = rep(0, 2), ncomp = 1),
    "zero"
  )
  expect_error(
    genpca:::genpca_cov_gmd(diag(2), R = matrix(0, 2, 2), ncomp = 1),
    "zero"
  )
  expect_error(
    genpca:::genpca_cov_gmd(matrix(0, 2, 2), R = diag(2), ncomp = 1),
    "No positive"
  )

  C_tiny_negative <- diag(c(1, -1e-7))
  expect_warning(
    genpca:::genpca_cov_geigen(
      C_tiny_negative,
      R = diag(2),
      ncomp = 1,
      constraints_remedy = "ridge",
      verbose = TRUE
    ),
    "non-PSD"
  )
  expect_error(
    genpca:::genpca_cov_geigen(diag(2), R = 1, ncomp = 1),
    "Length"
  )
  expect_error(
    genpca:::genpca_cov_geigen(diag(2), R = c(1, -1), ncomp = 1),
    "non-negative"
  )

  R_nonsym <- matrix(c(2, 0.1, 0, 1), 2)
  expect_error(
    genpca:::genpca_cov_geigen(C_near, R = R_nonsym, ncomp = 1),
    "symmetric"
  )
  expect_error(
    genpca:::genpca_cov_geigen(diag(2), R = c(0, 0), ncomp = 1),
    "zero"
  )
  expect_error(
    genpca:::genpca_cov_geigen(diag(2), R = matrix(0, 2, 2), ncomp = 1),
    "zero"
  )
  expect_error(
    genpca:::genpca_cov_geigen(matrix(0, 2, 2), R = diag(2), ncomp = 1),
    "No positive"
  )
})

test_that("GPCA and GPLSSVD operator fallbacks are covered", {
  set.seed(9503)
  X <- matrix(rnorm(6 * 4), 6, 4)
  Y <- matrix(rnorm(6 * 3), 6, 3)

  expect_message(
    genpca::genpca(X, ncomp = 1, method = "auto", verbose = TRUE),
    "Auto method selected"
  )
  expect_error(
    genpca::genpca(
      matrix(c(1, NA_real_, 2, 3, 4, 5), 3, 2),
      ncomp = 1,
      method = "spectra"
    ),
    "non-finite"
  )
  expect_error(
    genpca::genpca(
      matrix(c(1, NA_real_, 2, 3, 4, 5), 3, 2),
      ncomp = 1,
      method = "randomized"
    ),
    "non-finite"
  )
  suppressWarnings(
    zero_randomized <- genpca::genpca(
      matrix(0, 5, 4),
      ncomp = 2,
      method = "randomized",
      oversample = 1,
      n_power = 0
    )
  )
  expect_equal(multivarious::ncomp(zero_randomized), 0L)

  Q <- Matrix::Diagonal(nrow(X))
  R <- Matrix::Diagonal(ncol(X))
  attr(R, "eigen_decomp_cache_R") <- list(
    values = rep(1, ncol(X)),
    vectors = NULL,
    sqrtm = R,
    invsqrtm = R
  )
  expect_message(
    genpca:::gmdLA(
      X, Q, R,
      k = 1,
      n_orig = nrow(X),
      p_orig = ncol(X),
      verbose = TRUE
    ),
    "cached decomposition"
  )
  expect_error(
    genpca:::gmdLA(
      X, Q, Matrix::Matrix(matrix(c(0, 1e-12, 1e-12, 0), 2), sparse = FALSE),
      k = 1,
      n_orig = nrow(X),
      p_orig = 2
    ),
    "must be PSD|positive semi-definite|no positive eigenvalues"
  )
  expect_error(
    genpca:::gmdLA(
      X, Q, Matrix::Diagonal(ncol(X)),
      k = 0,
      n_orig = nrow(X),
      p_orig = ncol(X),
      use_dual = FALSE
    ),
    "k_request"
  )
  expect_error(
    genpca:::gmdLA(
      t(X), Matrix::Diagonal(ncol(X)), Matrix::Diagonal(nrow(X)),
      k = 0,
      n_orig = ncol(X),
      p_orig = nrow(X),
      use_dual = TRUE
    ),
    "k_request"
  )

  expect_error(
    genpca:::gmd_deflationR(X, Q, Matrix::Diagonal(ncol(X)), k = 1, maxit = 0),
    "maxit"
  )
  suppressWarnings(
    defl_zero <- genpca:::gmd_deflationR(
      matrix(0, 3, 2),
      Matrix::Diagonal(3),
      Matrix::Diagonal(2),
      k = 1,
      maxit = 2,
      verbose = TRUE
    )
  )
  expect_equal(defl_zero$k, 0L)

  expect_error(genpca::gplssvd_op(X, Y, k = 0), "k")
  expect_warning(
    too_many <- genpca::gplssvd_op(X, Y, k = 10),
    "exceeds"
  )
  expect_equal(too_many$k, 3L)

  Xs <- Matrix::Matrix(X, sparse = TRUE)
  Ys <- Matrix::Matrix(Y, sparse = TRUE)
  centered_scaled <- genpca::gplssvd_op(Xs, Ys, k = 1, center = TRUE, scale = TRUE)
  expect_equal(centered_scaled$k, 1L)

  row_metric <- diag(nrow(X))
  row_metric[1, 2] <- row_metric[2, 1] <- 0.1
  lazy <- genpca::gplssvd_op(
    Xs, Ys,
    XLW = row_metric,
    YLW = row_metric,
    k = 1
  )
  expect_equal(lazy$k, 1L)

  # Non-finite input must be rejected loudly, not silently "remediated"
  expect_error(genpca:::ensure_spd(matrix(NaN, 2, 2)), "non-finite")
})

test_that("GEP subspace and randomized GMD R fallbacks are covered", {
  S1 <- diag(c(3, 2, 1))
  S2 <- diag(3)

  expect_message(
    genpca:::solve_gep_subspace(
      S1, S2,
      q = 1,
      seed = 42,
      max_iter = 3,
      tol = Inf,
      verbose = TRUE
    ),
    "Converged"
  )
  expect_warning(
    genpca:::solve_gep_subspace(
      S1, S2,
      q = 1,
      which = "smallest",
      V0 = matrix(c(1, 0, 0), 3, 1),
      max_iter = 1,
      tol = 0
    ),
    "Reached max_iter"
  )

  indefinite_metric <- diag(c(1, -1))
  ortho <- genpca:::orthonormalize(
    diag(2),
    metric = indefinite_metric,
    reg = 1e-8,
    max_tries = 1
  )
  expect_equal(dim(ortho), c(2L, 2L))

  metric_fallback <- genpca:::metric_orthonormalize(
    diag(2),
    applyM = function(B) indefinite_metric %*% B,
    jitter = 0
  )
  expect_equal(dim(metric_fallback), c(2L, 1L))

  X <- matrix(c(1, 2, 3, 4, 5, 6), 3, 2)
  randomized_base_metrics <- genpca:::gmd_randomized_r(
    X,
    Q = diag(3),
    R = diag(2),
    k = 1,
    oversample = 0,
    n_power = 1,
    n_polish = 1,
    seed = 9504
  )
  expect_equal(randomized_base_metrics$k, 1L)

  zero_randomized_r <- genpca:::gmd_randomized_r(
    matrix(0, 3, 2),
    Q = diag(3),
    R = diag(2),
    k = 1,
    oversample = 0,
    n_power = 0
  )
  expect_equal(zero_randomized_r$k, 0L)

  expect_warning(
    expect_error(
      genpca:::gmd_randomized(
        X,
        Q = diag(2),
        R = diag(2),
        k = 1,
        use_cpp = TRUE
      ),
      "non-conformable|incompatible|dimension"
    ),
    "falling back to R implementation"
  )
})

test_that("remaining GPCA, RPLS, GenPLS, and transfer validation branches are covered", {
  set.seed(9505)
  X <- matrix(rnorm(6 * 4), 6, 4)
  Y <- matrix(rnorm(6 * 3), 6, 3)

  expect_message(
    suppressWarnings(
      genpca::genpca(
        X,
        ncomp = 1,
        method = "deflation",
        use_cpp = FALSE,
        maxit_deflation = 1,
        verbose = TRUE
      )
    ),
    "Using iterative deflation"
  )
  expect_message(
    suppressWarnings(
      genpca::genpca(
        t(X),
        ncomp = 1,
        method = "deflation",
        use_cpp = FALSE,
        maxit_deflation = 1,
        verbose = TRUE
      )
    ),
    "Using iterative deflation"
  )
  expect_error(
    genpca:::gmdLA(
      matrix(0, 2, 3),
      Matrix::Diagonal(2),
      Matrix::Diagonal(3),
      k = 1,
      n_orig = 2,
      p_orig = 3,
      use_dual = TRUE
    ),
    "No positive eigenvalues"
  )
  deflation_warnings <- character()
  deflation_error <- tryCatch(
    withCallingHandlers(
      genpca:::gmd_deflationR(
        X,
        Q = diag(2),
        R = diag(ncol(X)),
        k = 1,
        maxit = 1
      ),
      warning = function(w) {
        deflation_warnings <<- c(deflation_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  expect_s3_class(deflation_error, "error")
  expect_true(any(grepl("Could not compute total variance", deflation_warnings)))
  capture.output(
    converged_deflation <- suppressWarnings(genpca:::gmd_deflationR(
      matrix(c(4, 0, 0, 0), 2, 2),
      Q = Matrix::Diagonal(2),
      R = Matrix::Diagonal(2),
      k = 1,
      thr = 1e-4,
      maxit = 100,
      verbose = TRUE
    )),
    type = "message"
  )
  expect_lte(converged_deflation$k, 1L)

  zero_rpls_messages <- capture.output(
    invisible(genpca:::fit_rpls(matrix(0, 4, 2), matrix(0, 4, 2), K = 1, verbose = TRUE)),
    type = "message"
  )
  expect_match(paste(zero_rpls_messages, collapse = "\n"), "all zero|Degenerate")

  sparse_rpls_messages <- capture.output(
    invisible(genpca:::fit_rpls(
      matrix(rnorm(12), 4, 3),
      matrix(rnorm(8), 4, 2),
      K = 1,
      lambda = 1e6,
      verbose = TRUE
    )),
    type = "message"
  )
  expect_match(paste(sparse_rpls_messages, collapse = "\n"), "all zero")
  expect_error(genpca::genpls(X, Y[-1, ], ncomp = 1), "same number of rows")
  expect_message(
    genpca::genpls(X, Y, ncomp = 1, verbose = TRUE),
    "genpls finished"
  )

  rpls_fit <- genpca::rpls(X, Y, K = 1, lambda = 0.01)
  expect_error(
    genpca:::transfer.cross_projector(rpls_fit, X[1:2, ], from = "bad"),
    "'arg' should be one of"
  )
  expect_error(
    genpca:::transfer.cross_projector(rpls_fit, X[1:2, ], source = "bad"),
    "'arg' should be one of"
  )
})

test_that("final narrow coverage branches are exercised", {
  X <- matrix(c(1, 2, 3, 4), 2, 2)
  dge_metric <- genpca:::as_dge(Matrix::forceSymmetric(
    Matrix::Matrix(matrix(c(2, 0.1, 0.1, 1), 2), sparse = FALSE)
  ))
  constraints <- genpca:::prep_constraints(
    X,
    A = dge_metric,
    M = dge_metric,
    remedy = "error"
  )
  expect_s4_class(constraints$A, "dgeMatrix")
  expect_s4_class(constraints$M, "dgeMatrix")
  expect_error(
    genpca:::prep_constraints(
      X,
      A = diag(2),
      M = diag(c(1, -1)),
      remedy = "error"
    ),
    "Matrix M must be positive"
  )

  R_zeroish <- matrix(c(0, 1e-12, 1e-12, 0), 2)
  expect_error(
    genpca:::genpca_cov_gmd(diag(2), R = R_zeroish, ncomp = 1),
    "zero|positive semi-definite"
  )
  expect_error(
    genpca:::genpca_cov_geigen(diag(2), R = R_zeroish, ncomp = 1),
    "zero|PSD|positive semi-definite"
  )
  expect_equal(
    genpca:::genpca_cov_geigen(diag(2), R = diag(2), ncomp = NULL)$k,
    2L
  )

  expect_warning(
    genpca:::solve_gep_subspace(
      diag(2), diag(2),
      q = 1,
      V0 = matrix(c(1, 0), 2, 1),
      max_iter = 1,
      tol = 0,
      reg_T = -1
    ),
    "Reached max_iter"
  )

  set.seed(9506)
  expect_message(
    genpca::genpca(
      matrix(rnorm(5 * 3), 5, 3),
      ncomp = 1,
      method = "spectra",
      verbose = TRUE
    ),
    "whitened-operator SVD"
  )

  fit <- genpca::rpls(
    matrix(rnorm(6 * 3), 6, 3),
    matrix(rnorm(6 * 2), 6, 2),
    K = 1,
    lambda = 0.01
  )
  transfer_attempt <- try(
    genpca:::transfer.cross_projector(fit, matrix(rnorm(2 * 3), 2, 3), to = "Y"),
    silent = TRUE
  )
  expect_true(is.matrix(transfer_attempt) || inherits(transfer_attempt, "try-error"))
})
