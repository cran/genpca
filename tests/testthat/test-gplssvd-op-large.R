# The toy tests in test-gplssvd_op.R all fall into the small-dense fallback
# (ncol <= 64). These tests exercise the iterative operator paths used for
# larger problems: materialized whitened blocks (dense and sparse-preserving),
# the lazy chain for sparse data with dense-general metrics, and both SVD
# backends -- checked against an explicit dense whitening reference.

suppressPackageStartupMessages(library(Matrix))

sym_sqrt <- function(W) {
  e <- eigen(as.matrix(W), symmetric = TRUE)
  e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors))
}

half_mat <- function(W, m) {
  if (is.null(W)) {
    diag(m)
  } else if (is.numeric(W) && is.null(dim(W))) {
    diag(sqrt(pmax(W, 0)), m)
  } else if (inherits(W, "diagonalMatrix")) {
    diag(sqrt(pmax(as.numeric(Matrix::diag(W)), 0)), m)
  } else {
    sym_sqrt(W)
  }
}

dense_ref <- function(X, Y, MX = NULL, MY = NULL, WX = NULL, WY = NULL, k) {
  Xe <- half_mat(MX, nrow(X)) %*% as.matrix(X) %*% half_mat(WX, ncol(X))
  Ye <- half_mat(MY, nrow(Y)) %*% as.matrix(Y) %*% half_mat(WY, ncol(Y))
  sv <- svd(crossprod(Xe, Ye))
  list(d = sv$d[seq_len(k)],
       u = sv$u[, seq_len(k), drop = FALSE],
       v = sv$v[, seq_len(k), drop = FALSE])
}

expect_matches_ref <- function(op, ref, tol = 1e-6) {
  testthat::expect_equal(op$d, ref$d, tolerance = tol)
  # compare subspaces columnwise via |<u_i, u_ref_i>| = 1 (sign-invariant)
  cu <- abs(colSums(as.matrix(op$u) * ref$u))
  cv <- abs(colSums(as.matrix(op$v) * ref$v))
  testthat::expect_true(all(abs(cu - 1) < tol))
  testthat::expect_true(all(abs(cv - 1) < tol))
}

set.seed(101)
N <- 150
I <- 80
J <- 70
k <- 3
X <- matrix(rnorm(N * I), N, I)
Y <- matrix(rnorm(N * J), N, J)
Xs <- rsparsematrix(N, I, 0.08)
Ys <- rsparsematrix(N, J, 0.08)
wr <- runif(N)
wcx <- runif(I)
wcy <- runif(J)
mk_spd <- function(m) {
  A <- matrix(rnorm(m * m, sd = 0.3), m, m)
  crossprod(A) + diag(m)
}
Mg <- mk_spd(N)
Mg2 <- mk_spd(N)
Wg <- mk_spd(I)

testthat::test_that("eigencore path: identity and diagonal metrics (materialized)", {
  expect_matches_ref(gplssvd_op(X, Y, k = k), dense_ref(X, Y, k = k))
  op <- gplssvd_op(X, Y, XLW = Diagonal(x = wr), YLW = Diagonal(x = wr),
                   XRW = wcx, YRW = wcy, k = k)
  ref <- dense_ref(X, Y, MX = Diagonal(x = wr), MY = Diagonal(x = wr),
                   WX = wcx, WY = wcy, k = k)
  expect_matches_ref(op, ref)
})

testthat::test_that("eigencore path: general SPD metrics on dense data", {
  op <- gplssvd_op(X, Y, XLW = Mg, YLW = Mg, XRW = Wg, k = k)
  expect_matches_ref(op, dense_ref(X, Y, MX = Mg, MY = Mg, WX = Wg, k = k))
  # distinct row metrics on each block
  op2 <- gplssvd_op(X, Y, XLW = Mg, YLW = Mg2, k = k)
  expect_matches_ref(op2, dense_ref(X, Y, MX = Mg, MY = Mg2, k = k))
})

testthat::test_that("eigencore path: sparse data stays on the lazy chain with general metrics", {
  # shared row metric (fused middle) and distinct row metrics
  op <- gplssvd_op(Xs, Ys, XLW = Mg, YLW = Mg, k = k)
  expect_matches_ref(op, dense_ref(Xs, Ys, MX = Mg, MY = Mg, k = k))
  op2 <- gplssvd_op(Xs, Ys, XLW = Mg, YLW = Mg2, k = k)
  expect_matches_ref(op2, dense_ref(Xs, Ys, MX = Mg, MY = Mg2, k = k))
})

testthat::test_that("eigencore path: sparse data with diagonal metrics (sparse materialization)", {
  op <- gplssvd_op(Xs, Ys, XLW = Diagonal(x = wr), YLW = Diagonal(x = wr),
                   XRW = wcx, YRW = wcy, k = k)
  ref <- dense_ref(Xs, Ys, MX = Diagonal(x = wr), MY = Diagonal(x = wr),
                   WX = wcx, WY = wcy, k = k)
  expect_matches_ref(op, ref)
})

testthat::test_that("irlba backend works beyond the dense fallback threshold", {
  testthat::skip_if_not_installed("irlba")
  op <- gplssvd_op(X, Y, k = k, svd_backend = "irlba")
  expect_matches_ref(op, dense_ref(X, Y, k = k))
  op2 <- gplssvd_op(Xs, Ys, XLW = Mg, YLW = Mg2, k = k, svd_backend = "irlba")
  expect_matches_ref(op2, dense_ref(Xs, Ys, MX = Mg, MY = Mg2, k = k))
})

testthat::test_that("derived quantities match explicit formulas on the iterative path", {
  op <- gplssvd_op(X, Y, XLW = Diagonal(x = wr), YLW = Diagonal(x = wr),
                   XRW = wcx, YRW = wcy, k = k)
  ref <- dense_ref(X, Y, MX = Diagonal(x = wr), MY = Diagonal(x = wr),
                   WX = wcx, WY = wcy, k = k)
  sgn <- sign(colSums(as.matrix(op$u) * ref$u))
  p_ref <- (1 / sqrt(wcx)) * (ref$u %*% diag(sgn, k))
  q_ref <- (1 / sqrt(wcy)) * (ref$v %*% diag(sgn, k))
  fi_ref <- wcx * p_ref %*% diag(ref$d, k)
  fj_ref <- wcy * q_ref %*% diag(ref$d, k)
  lx_ref <- sqrt(wr) * (X %*% (wcx * p_ref))
  ly_ref <- sqrt(wr) * (Y %*% (wcy * q_ref))
  testthat::expect_equal(as.matrix(op$p), p_ref, tolerance = 1e-6, ignore_attr = TRUE)
  testthat::expect_equal(as.matrix(op$q), q_ref, tolerance = 1e-6, ignore_attr = TRUE)
  testthat::expect_equal(as.matrix(op$fi), fi_ref, tolerance = 1e-6, ignore_attr = TRUE)
  testthat::expect_equal(as.matrix(op$fj), fj_ref, tolerance = 1e-6, ignore_attr = TRUE)
  testthat::expect_equal(as.matrix(op$lx), lx_ref, tolerance = 1e-6, ignore_attr = TRUE)
  testthat::expect_equal(as.matrix(op$ly), ly_ref, tolerance = 1e-6, ignore_attr = TRUE)
})

testthat::test_that("genpls keeps sparse inputs sparse with pass() and matches dense fit", {
  testthat::skip_if_not_installed("multivarious")
  Mw <- Diagonal(x = wr)
  fit_sp <- genpls(Xs, Ys, Mx = Mw, My = Mw, ncomp = k)
  fit_de <- genpls(as.matrix(Xs), as.matrix(Ys), Mx = Mw, My = Mw, ncomp = k)
  testthat::expect_equal(fit_sp$d, fit_de$d, tolerance = 1e-6)
  sgn <- sign(colSums(fit_sp$vx * fit_de$vx))
  testthat::expect_lt(max(abs(fit_sp$vx - fit_de$vx %*% diag(sgn, k))), 1e-6)
  testthat::expect_lt(max(abs(fit_sp$vy - fit_de$vy %*% diag(sgn, k))), 1e-6)

  # the sparse-fitted projector must project new data identically
  set.seed(142)
  newX <- matrix(rnorm(5 * I), 5, I)
  pr_sp <- multivarious::project(fit_sp, newX, source = "X")
  pr_de <- multivarious::project(fit_de, newX, source = "X")
  testthat::expect_lt(max(abs(as.matrix(pr_sp) - as.matrix(pr_de) %*% diag(sgn, k))), 1e-6)
})

testthat::test_that("genpls with centering preprocessor densifies sparse input gracefully", {
  testthat::skip_if_not_installed("multivarious")
  fit_c <- genpls(Xs, Ys, ncomp = k,
                  preproc_x = multivarious::center(),
                  preproc_y = multivarious::center())
  fit_cd <- genpls(as.matrix(Xs), as.matrix(Ys), ncomp = k,
                   preproc_x = multivarious::center(),
                   preproc_y = multivarious::center())
  testthat::expect_equal(fit_c$d, fit_cd$d, tolerance = 1e-6)
})

testthat::test_that("mismatched metric dimensions error clearly", {
  testthat::expect_error(gplssvd_op(X, Y, XLW = runif(N - 1), k = k),
                         "length")
  testthat::expect_error(gplssvd_op(X, Y, XRW = mk_spd(I + 2), k = k),
                         "dimension")
})
