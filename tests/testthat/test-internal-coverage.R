test_that("partial eigen and metric sqrt helpers cover diagonal and dense paths", {
  M_diag <- Matrix::Diagonal(3, x = c(4, 1, 0.25))
  expect_true(genpca:::is_identity_or_diag(M_diag))
  expect_false(genpca:::is_identity_or_diag(diag(3)))

  eig <- genpca:::partial_eig_once(as.matrix(M_diag), k = 2)
  expect_equal(eig$lam, c(4, 1), tolerance = 1e-12)
  expect_equal(dim(eig$Q), c(3L, 2L))

  M_big <- diag(as.numeric(seq(60, 1)))
  eig_big <- genpca:::partial_eig_once(M_big, k = 2, tol = 1e-8)
  expect_equal(sort(eig_big$lam, decreasing = TRUE), c(60, 59), tolerance = 1e-8)

  adapt <- genpca:::decide_adaptive_rank(as.matrix(M_diag), var_threshold = 0.90, max_k = 3)
  expect_equal(length(adapt$lam), 2L)

  approx_fixed <- genpca:::partial_eig_approx(as.matrix(M_diag), user_rank = 2)
  approx_adapt <- genpca:::partial_eig_approx(as.matrix(M_diag), user_rank = NA, max_k = 3)
  expect_equal(length(approx_fixed$lams), 2L)
  expect_true(length(approx_adapt$lams) >= 1L)

  X <- matrix(c(1, 2, 3, 4, 5, 6), 3, 2)
  expect_equal(genpca:::build_sqrt_mult(NULL, NULL)(X), X)
  expect_equal(as.matrix(genpca:::build_sqrt_mult(M_diag, NULL)(X)),
               as.matrix(Matrix::Diagonal(3, x = c(2, 1, 0.5)) %*% X))
  expect_equal(as.matrix(genpca:::build_invsqrt_mult(M_diag, NULL)(X)),
               as.matrix(Matrix::Diagonal(3, x = c(0.5, 1, 2)) %*% X))

  M_dense <- matrix(c(2, 1, 1, 2), 2, 2)
  X2 <- matrix(c(1, 3, 2, 4), 2, 2)
  ed <- eigen(M_dense, symmetric = TRUE)
  S <- ed$vectors %*% (sqrt(ed$values) * t(ed$vectors))
  IS <- ed$vectors %*% ((1 / sqrt(ed$values)) * t(ed$vectors))
  expect_equal(as.matrix(genpca:::build_sqrt_mult(M_dense, 2)(X2)),
               S %*% X2, tolerance = 1e-10)
  expect_equal(as.matrix(genpca:::build_invsqrt_mult(M_dense, 2)(X2)),
               IS %*% X2, tolerance = 1e-10)

  M_rank1 <- matrix(c(1, 1, 1, 1), 2, 2)
  f_rank1 <- genpca:::build_invsqrt_mult(M_rank1, 2, eps = 1e-12)
  out_rank1 <- f_rank1(cbind(c(1, 1), c(1, -1)))
  expect_equal(out_rank1[, 2], c(0, 0), tolerance = 1e-10)

  expect_equal(genpca:::row_transform(X2, function(z) z + 1), X2 + 1)
  expect_equal(genpca:::col_transform(X2, function(z) 2 * z), 2 * X2)
})

test_that("weight operators apply identity, diagonal, and dense transformations", {
  x <- matrix(c(1, 2, 3, 4), 2, 2)
  expect_equal(genpca:::as_weight_operator(NULL)(x), x)

  Wd <- Matrix::Diagonal(2, x = c(4, 9))
  expect_equal(as.matrix(genpca:::as_weight_operator(Wd)(x)),
               as.matrix(Matrix::Diagonal(2, x = c(4, 9)) %*% x))
  expect_equal(as.matrix(genpca:::as_weight_operator(Wd, sqrt = TRUE)(x)),
               as.matrix(Matrix::Diagonal(2, x = c(2, 3)) %*% x))
  expect_equal(as.matrix(genpca:::as_weight_operator(Wd, inverse = TRUE)(x)),
               as.matrix(Matrix::Diagonal(2, x = c(1 / 4, 1 / 9)) %*% x))
  expect_equal(as.matrix(genpca:::as_weight_operator(Wd, sqrt = TRUE, inverse = TRUE)(x)),
               as.matrix(Matrix::Diagonal(2, x = c(1 / 2, 1 / 3)) %*% x))

  W <- matrix(c(2, 1, 1, 2), 2, 2)
  ed <- eigen(W, symmetric = TRUE)
  S <- ed$vectors %*% (sqrt(ed$values) * t(ed$vectors))
  IS <- ed$vectors %*% ((1 / sqrt(ed$values)) * t(ed$vectors))
  expect_equal(as.matrix(genpca:::as_weight_operator(W, transpose = TRUE)(x)),
               W %*% x, tolerance = 1e-10)
  expect_equal(as.matrix(genpca:::as_weight_operator(W, sqrt = TRUE)(x)),
               S %*% x, tolerance = 1e-10)
  expect_equal(as.matrix(genpca:::as_weight_operator(W, sqrt = TRUE, inverse = TRUE)(x)),
               IS %*% x, tolerance = 1e-10)
  expect_equal(as.matrix(genpca:::as_weight_operator(W, inverse = TRUE)(x)),
               solve(W, x), tolerance = 1e-10)
})

test_that("constraint utilities coerce matrices and repair indefinite constraints", {
  expect_true(genpca:::is_spd(diag(2)))
  expect_false(genpca:::is_spd(matrix(c(1, 2, 0, 1), 2, 2)))
  expect_false(genpca:::is_spd(matrix(c(1, 0, 0, -1), 2, 2)))

  dg <- Matrix::Diagonal(3, x = c(1, 2, 3))
  dgc <- genpca:::as_dgc(dg)
  expect_s4_class(dgc, "dgCMatrix")
  expect_identical(genpca:::as_dgc(dgc), dgc)

  dense <- Matrix::Matrix(matrix(c(2, 1, 1, 2), 2, 2), sparse = FALSE)
  dge <- genpca:::as_dge(dense)
  expect_s4_class(dge, "dgeMatrix")
  expect_identical(genpca:::as_dge(dge), dge)

  M <- matrix(c(0, 2, 2, 0), 2, 2)
  M_spd <- genpca:::ensure_spd(M, tol = 1e-6)
  expect_true(genpca:::is_spd(M_spd))
  expect_true(Matrix::isSymmetric(M_spd))
})

test_that("GMD cache digests, reuses factors, evicts least-recent entries, and clears", {
  genpca::gmd_clear_cache()
  expect_identical(genpca::gmd_clear_cache(), TRUE)

  A <- Matrix::Matrix(crossprod(matrix(c(1, 2, 3, 4), 2, 2)) + diag(2), sparse = FALSE)
  digest_dense <- genpca:::.digest_dense_matrix(A)
  digest_sparse <- genpca:::.digest_sparse_matrix(Matrix::Matrix(A, sparse = TRUE))
  expect_type(digest_dense, "character")
  expect_type(digest_sparse, "character")
  expect_error(genpca:::.digest_sparse_matrix(A))

  L1 <- genpca:::get_chol_lower_dense(A)
  L2 <- genpca:::get_chol_lower_dense(A)
  expect_equal(L1, L2)
  expect_true(length(ls(envir = genpca:::.gmd_cache, all.names = TRUE)) >= 1L)

  for (i in seq_len(18)) {
    B <- diag(c(i + 1, i + 2))
    genpca:::get_chol_lower_dense(B)
  }
  expect_lte(length(ls(envir = genpca:::.gmd_cache, all.names = TRUE)), 16L)
  expect_type(genpca:::.evict_lru(), "character")

  genpca::gmd_clear_cache()
  expect_length(ls(envir = genpca:::.gmd_cache, all.names = TRUE), 0L)
  expect_length(ls(envir = genpca:::.gmd_cache_times, all.names = TRUE), 0L)
  expect_null(genpca:::.evict_lru())
})

test_that("GMD fast R helpers validate options and cover zero-rank and randomized fallback paths", {
  expect_false(genpca:::should_use_topk(1, 100, topk = FALSE, auto_topk = TRUE, topk_ratio = 0.1, topk_min_dim = 20))
  expect_true(genpca:::should_use_topk(1, 10, topk = TRUE, auto_topk = FALSE, topk_ratio = 0.1, topk_min_dim = 20))
  expect_false(genpca:::should_use_topk(1, 10, topk = TRUE, auto_topk = TRUE, topk_ratio = 0.1, topk_min_dim = 20))
  expect_true(genpca:::should_use_topk(5, 100, topk = TRUE, auto_topk = TRUE, topk_ratio = 0.1, topk_min_dim = 20))

  zero_diag <- genpca:::gmd_fast_cpp(
    matrix(0, 4, 3),
    Matrix::Diagonal(4),
    Matrix::Diagonal(3),
    k = 2,
    tol = 1e-9
  )
  expect_equal(zero_diag$k, 0L)
  expect_equal(dim(zero_diag$u), c(4L, 0L))

  set.seed(126)
  X <- matrix(rnorm(30), 6, 5)
  Q <- Matrix::Diagonal(6)
  R <- Matrix::Diagonal(5)
  applyI <- function(z) z
  expect_equal(dim(genpca:::metric_orthonormalize(matrix(numeric(0), 4, 0), applyI)), c(4L, 0L))
  B <- genpca:::metric_orthonormalize(X[, 1:2], applyI)
  expect_equal(crossprod(B), diag(2), tolerance = 1e-8)

  signs1 <- genpca:::random_sign_matrix(4, 3, seed = 10)
  signs2 <- genpca:::random_sign_matrix(4, 3, seed = 10)
  expect_equal(signs1, signs2)
  expect_true(all(signs1 %in% c(-1, 1)))

  expect_error(genpca:::gmd_randomized_r(X, Q, R, k = 0), "positive integer")
  expect_error(genpca:::gmd_randomized_r(X, Q, R, k = 2, oversample = -1), "non-negative integer")
  expect_error(genpca:::gmd_randomized_r(X, Q, R, k = 2, n_power = -1), "non-negative integer")
  expect_error(genpca:::gmd_randomized_r(X, Q, R, k = 2, n_polish = -1), "non-negative integer")
  expect_error(genpca:::gmd_randomized_r(X, Q, R, k = 2, polish_tol = -1), "non-negative")

  rand <- genpca:::gmd_randomized(
    X, Q, R, k = 2,
    oversample = 2,
    n_power = 1,
    n_polish = 1,
    polish_tol = Inf,
    seed = 99,
    use_cpp = FALSE
  )
  expect_equal(rand$k, 2L)
  expect_equal(as.matrix(crossprod(rand$u)), diag(2), tolerance = 1e-6)
  expect_equal(as.matrix(crossprod(rand$v)), diag(2), tolerance = 1e-6)

  rand_zero <- genpca:::gmd_randomized_r(matrix(0, 4, 3), Matrix::Diagonal(4), Matrix::Diagonal(3), k = 2)
  expect_equal(rand_zero$k, 0L)

  pol0 <- genpca:::gmd_randomized_polish(X, applyI, applyI, rand$u, rand$v, iters = 0)
  expect_equal(pol0$d, numeric(2))
  pol1 <- genpca:::gmd_randomized_polish(X, applyI, applyI, rand$u, rand$v, iters = 2, polish_tol = Inf)
  expect_equal(length(pol1$d), 2L)
})

test_that("GMD fast dispatch covers dense, sparse, cached, uncached, primal, and dual paths", {
  set.seed(123)
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

  X_primal <- matrix(rnorm(8 * 5), 8, 5)
  Q_sp <- make_sparse_spd(8)
  R_dn <- Matrix::Matrix(make_spd(5), sparse = FALSE)
  res_primal_sp <- genpca:::gmd_fast_cpp(X_primal, Q_sp, R_dn, k = 2, diag_fast = FALSE)
  expect_equal(res_primal_sp$k, 2L)

  res_uncached <- genpca:::gmd_fast_cpp(
    X_primal,
    Matrix::Matrix(make_spd(8), sparse = FALSE),
    Matrix::Matrix(make_spd(5), sparse = FALSE),
    k = 2,
    cache = FALSE,
    diag_fast = FALSE,
    topk = FALSE
  )
  expect_equal(res_uncached$k, 2L)

  X_dual <- matrix(rnorm(5 * 8), 5, 8)
  Q_dn <- Matrix::Matrix(make_spd(5), sparse = FALSE)
  R_sp <- make_sparse_spd(8)
  res_dual_sp <- genpca:::gmd_fast_cpp(X_dual, Q_dn, R_sp, k = 2, diag_fast = FALSE)
  expect_equal(res_dual_sp$k, 2L)

  res_dual_uncached <- genpca:::gmd_fast_cpp(
    X_dual,
    Matrix::Matrix(make_spd(5), sparse = FALSE),
    Matrix::Matrix(make_spd(8), sparse = FALSE),
    k = 2,
    cache = FALSE,
    diag_fast = FALSE,
    topk = FALSE
  )
  expect_equal(res_dual_uncached$k, 2L)

  expect_error(genpca:::gmd_fast_cpp(X_primal, Matrix::Diagonal(8), Matrix::Diagonal(5), k = 0), "k must")
  expect_warning(
    clipped <- genpca:::gmd_fast_cpp(X_primal, Matrix::Diagonal(8), Matrix::Diagonal(5), k = 99),
    "exceeds"
  )
  expect_lte(clipped$k, 5L)
  expect_error(genpca:::gmd_fast_cpp(X_primal, Matrix::Diagonal(8), Matrix::Diagonal(5), k = 2, maxit = 0), "maxit")
  expect_error(genpca:::gmd_fast_cpp(X_primal, Matrix::Diagonal(8), Matrix::Diagonal(5), k = 2, topk_ratio = 0), "topk_ratio")
  expect_error(genpca:::gmd_fast_cpp(X_primal, Matrix::Diagonal(8), Matrix::Diagonal(5), k = 2, topk_min_dim = 1), "topk_min_dim")
})

test_that("randomized C++ entry points cover zero and mixed sparse/dense metric combinations", {
  set.seed(456)
  X <- matrix(rnorm(7 * 5), 7, 5)
  Qd <- diag(7)
  Rd <- diag(5)
  Qs <- Matrix::Diagonal(7)
  Rs <- Matrix::Diagonal(5)

  zero <- genpca:::gmd_randomized_cpp_dn(X, Qd, Rd, k = 0)
  expect_equal(zero$k, 0L)

  dense <- genpca:::gmd_randomized_cpp_dn(X, Qd, Rd, k = 2, oversample = 1, n_power = 1, n_polish = 1, seed = 7)
  both_sparse <- genpca:::gmd_randomized_cpp_sp(X, Qs, Rs, k = 2, oversample = 1, n_power = 1, n_polish = 1, seed = 7)
  q_sparse <- genpca:::gmd_randomized_cpp_qsp_rdn(X, Qs, Rd, k = 2, oversample = 1, n_power = 1, seed = 7)
  r_sparse <- genpca:::gmd_randomized_cpp_qdn_rsp(X, Qd, Rs, k = 2, oversample = 1, n_power = 1, seed = 7)

  expect_equal(dense$k, 2L)
  expect_equal(both_sparse$k, 2L)
  expect_equal(q_sparse$k, 2L)
  expect_equal(r_sparse$k, 2L)
  expect_equal(dense$d, both_sparse$d, tolerance = 1e-6)
})

test_that("GEP low-level helpers handle successful and failed regularization paths", {
  S <- Matrix::Diagonal(3, x = c(1, 2, 3))
  fact <- genpca:::factor_mat(S, reg = 1e-8)
  expect_named(fact, c("ch", "final_reg"))
  expect_error(genpca:::factor_mat(-Matrix::Diagonal(2), reg = 1e-12, max_tries = 1),
               "Unable to factor")

  X0 <- matrix(numeric(0), 3, 0)
  expect_equal(dim(genpca:::orthonormalize(X0)), c(3L, 0L))

  set.seed(266)
  X <- matrix(rnorm(9), 3, 3)
  Q <- qr.Q(qr(X))
  expect_equal(crossprod(genpca:::orthonormalize(X)), diag(3), tolerance = 1e-8)

  M <- Matrix::Diagonal(3, x = c(1, 2, 3))
  B <- genpca:::orthonormalize(Q[, 1:2], metric = M)
  expect_equal(as.matrix(crossprod(B, M %*% B)), diag(2), tolerance = 1e-6)

  V0 <- diag(3)[, 1:2]
  sol <- genpca:::solve_gep_subspace(
    Matrix::Diagonal(3, x = c(5, 3, 1)),
    Matrix::Diagonal(3),
    q = 2,
    which = "largest",
    max_iter = 20,
    V0 = V0,
    tol = 1e-8
  )
  expect_equal(sort(sol$values, decreasing = TRUE), c(5, 3), tolerance = 1e-4)
})

test_that("rpls internal fitter covers validation, degenerate, verbose, and ridge-GPLS edge paths", {
  set.seed(1201)
  X <- matrix(rnorm(8 * 4), 8, 4)
  Y <- matrix(rnorm(8 * 3), 8, 3)

  expect_error(genpca:::fit_rpls(X, Y, K = 3, lambda = c(0.1, 0.2)), "length >= K")
  expect_error(genpca:::fit_rpls(X, Y, K = 1, Q = diag(3)), "Q must be p x p")

  fit_vec_lambda <- genpca:::fit_rpls(
    X, Y,
    K = 2,
    lambda = c(0.05, 0.1),
    penalty = "l1",
    maxiter = 20,
    verbose = TRUE
  )
  expect_equal(fit_vec_lambda$num_components, 2L)

  fit_degenerate <- genpca:::fit_rpls(
    matrix(0, 8, 4),
    Y,
    K = 2,
    lambda = 0.1,
    penalty = "l1",
    verbose = TRUE
  )
  expect_equal(fit_degenerate$num_components, 0L)
  expect_equal(dim(fit_degenerate$V), c(4L, 0L))

  fit_zero_loading <- genpca:::fit_rpls(
    X, Y,
    K = 2,
    lambda = 1e6,
    penalty = "l1",
    verbose = TRUE
  )
  expect_equal(fit_zero_loading$num_components, 0L)

  expect_error(
    expect_warning(
      genpca:::fit_rpls(
        X, Y,
        K = 1,
        lambda = 0.1,
        penalty = "ridge",
        Q = -diag(4)
      ),
      "Cholesky decomposition failed"
    ),
    "Cannot proceed"
  )

  expect_warning(
    fit_ridge_nonneg <- genpca:::fit_rpls(
      X, Y,
      K = 1,
      lambda = 0.1,
      penalty = "ridge",
      Q = diag(4),
      nonneg = TRUE
    ),
    "nonneg option is ignored"
  )
  expect_equal(fit_ridge_nonneg$num_components, 1L)
})

test_that("remaining GMD wrappers cover rank filtering, zero-size randomized input, and dense C++ randomized dispatch", {
  set.seed(1202)
  X_rank1 <- tcrossprod(rnorm(6), rnorm(4))
  filtered <- genpca:::gmd_fast_cpp(
    X_rank1,
    Matrix::Matrix(diag(6), sparse = FALSE),
    Matrix::Matrix(diag(4), sparse = FALSE),
    k = 4,
    tol = 1e-6,
    cache = FALSE,
    diag_fast = FALSE,
    topk = FALSE
  )
  expect_equal(filtered$k, 1L)

  empty_rand <- genpca:::gmd_randomized_r(
    matrix(numeric(0), 0, 3),
    Matrix::Matrix(0, 0, 0),
    Matrix::Diagonal(3),
    k = 1
  )
  expect_equal(empty_rand$k, 0L)

  X <- matrix(rnorm(7 * 5), 7, 5)
  dense_rand <- genpca:::gmd_randomized(
    X,
    Matrix::Matrix(diag(7), sparse = FALSE),
    Matrix::Matrix(diag(5), sparse = FALSE),
    k = 2,
    oversample = 1,
    n_power = 1,
    n_polish = 1,
    polish_tol = Inf,
    seed = 1202,
    use_cpp = TRUE
  )
  expect_equal(dense_rand$k, 2L)

  expect_s4_class(genpca:::as_dgc(matrix(c(1, 0, 0, 1), 2, 2)), "dgCMatrix")
})
