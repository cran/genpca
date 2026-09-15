# These tests compare the memory-safe operator SVD (gplssvd_op)
# against the explicit "full SVD on the whitened cross-product" for toy cases.

testthat::test_that("setup: Matrix available", {
  testthat::skip_if_not_installed("Matrix")
})

suppressPackageStartupMessages({
  library(Matrix)
})

# ---- Helpers ---------------------------------------------------------------

set.seed(1)

# Build small SPD (positive definite) matrix; optionally return sparse.
make_spd <- function(n, dense = TRUE) {
  A <- matrix(rnorm(n * n, sd = 0.6), n, n)
  S <- crossprod(A) + diag(n) * 0.2
  if (dense) Matrix(S, sparse = FALSE) else as(forceSymmetric(Matrix(S, sparse = TRUE)), "dsCMatrix")
}

# Center/scale to match the algorithm in gplssvd_op()
center_scale_like_impl <- function(A, do_center, do_scale) {
  A <- as.matrix(A)
  if (is.null(dim(A))) A <- matrix(A, nrow = length(A), ncol = 1)
  cen <- if (isTRUE(do_center)) colMeans(A) else rep(0, ncol(A))
  if (isTRUE(do_center)) {
    A <- A - matrix(rep(cen, each = nrow(A)), nrow(A), byrow = FALSE)
  }
  if (isTRUE(do_scale)) {
    s <- sqrt(base::colSums(A^2) / pmax(nrow(A) - 1, 1))
    s[s == 0] <- 1
    A <- A %*% diag(1 / as.numeric(s), ncol(A))
  } else {
    s <- rep(1, ncol(A))
  }
  list(A = Matrix::Matrix(A, sparse = FALSE), center = cen, scale = s)
}

# For diagonal metric represented by numeric vector or diagonalMatrix,
# return explicit sqrt and inverse sqrt matrices.
make_sqrt_mats <- function(W, n) {
  if (is.null(W)) {
    list(sqrt = Diagonal(n), invsqrt = Diagonal(n), full = Diagonal(n))
  } else if (is.numeric(W) && length(W) == n) {
    d <- pmax(as.numeric(W), 0)
    ds <- sqrt(d)
    list(sqrt = Diagonal(x = ds),
         invsqrt = Diagonal(x = ifelse(ds > 0, 1 / ds, 0)),
         full = Diagonal(x = d))
  } else if (inherits(W, "diagonalMatrix")) {
    d <- pmax(as.numeric(diag(W)), 0)
    ds <- sqrt(d)
    list(sqrt = Diagonal(x = ds),
         invsqrt = Diagonal(x = ifelse(ds > 0, 1 / ds, 0)),
         full = W)
  } else {
    # General PSD: symmetric sqrt via eigen
    Wd <- as.matrix(forceSymmetric(Matrix(W)))
    es <- eigen(Wd, symmetric = TRUE)
    Q  <- es$vectors
    lam <- pmax(es$values, 0)
    S  <- Q %*% diag(sqrt(lam), nrow = length(lam)) %*% t(Q)
    IS <- Q %*% diag(ifelse(lam > 0, 1 / sqrt(lam), 0), nrow = length(lam)) %*% t(Q)
    list(
      sqrt    = Matrix(S,  sparse = FALSE),
      invsqrt = Matrix(IS, sparse = FALSE),
      full    = Matrix(Wd, sparse = FALSE)
    )
  }
}

# Dense reference implementation of GPLSSVD via explicit whitening (eqs. 11–14 in Beaton 2020).
dense_gplssvd_ref <- function(X, Y,
                              MX = NULL, MY = NULL,
                              WX = NULL, WY = NULL,
                              k = NULL,
                              center = FALSE, scale = FALSE) {
  N <- nrow(X)
  I <- ncol(X)
  J <- ncol(Y)
  stopifnot(nrow(Y) == N)

  # center/scale as in the operator impl
  Xcs <- center_scale_like_impl(X, center, scale)
  X <- Xcs$A
  Ycs <- center_scale_like_impl(Y, center, scale)
  Y <- Ycs$A

  # sqrt and invsqrt metric matrices (explicit)
  MXm <- make_sqrt_mats(MX, N)
  MYm <- make_sqrt_mats(MY, N)
  WXm <- make_sqrt_mats(WX, I)
  WYm <- make_sqrt_mats(WY, J)

  # Whitened matrices
  Xe <- MXm$sqrt %*% X %*% WXm$sqrt
  Ye <- MYm$sqrt %*% Y %*% WYm$sqrt

  # Core S = t(Xe) %*% Ye (small I x J)
  S <- crossprod(Xe, Ye)

  # SVD of S
  sv <- svd(as.matrix(S))
  if (!is.null(k)) {
    k <- min(k, length(sv$d))
    sv$u <- sv$u[, seq_len(k), drop = FALSE]
    sv$v <- sv$v[, seq_len(k), drop = FALSE]
    sv$d <- sv$d[seq_len(k)]
  }

  # Generalized singular vectors and scores (eq. 12 and 14)
  p <- WXm$invsqrt %*% sv$u
  q <- WYm$invsqrt %*% sv$v

  Fi <- WXm$full %*% (p %*% diag(sv$d, nrow = length(sv$d)))
  Fj <- WYm$full %*% (q %*% diag(sv$d, nrow = length(sv$d)))

  # Latent variables
  Lx <- MXm$sqrt %*% (X %*% (WXm$full %*% p))
  Ly <- MYm$sqrt %*% (Y %*% (WYm$full %*% q))

  list(d = sv$d, u = sv$u, v = sv$v, p = p, q = q, fi = Fi, fj = Fj, lx = Lx, ly = Ly,
       S = S,
       center = list(X = Xcs$center, Y = Ycs$center),
       scale  = list(X = Xcs$scale,  Y = Ycs$scale))
}

# align_signs() comes from helper-align.R

# Quick numeric check helper
expect_mats_equal <- function(A, B, tol = 1e-7) {
  testthat::expect_true(all(dim(A) == dim(B)))
  testthat::expect_lt(norm(as.matrix(A - B), "F"), tol * (1 + norm(as.matrix(A), "F")))
}

# Wrapper to run operator version with a chosen backend
run_op <- function(X, Y, XLW, YLW, XRW, YRW, k, center, scale, backend) {
  gplssvd_op(X, Y, XLW = XLW, YLW = YLW, XRW = XRW, YRW = YRW,
             k = k, center = center, scale = scale,
             svd_backend = backend,
             svd_opts = list(tol = 1e-9, maxitr = 5000))
}

# ---- Toy data --------------------------------------------------------------

N <- 12
I <- 5
J <- 4
X0 <- matrix(rnorm(N * I), N, I)
Y0 <- matrix(rnorm(N * J), N, J)

# ---- 1) Identity constraints (canonical PLS / PLS-SVD) --------------------

for (backend in c("eigencore", "irlba")) {
  testthat::test_that(paste("Identity metrics match full SVD (backend =", backend, ")"), {
    if (backend == "irlba") testthat::skip_if_not_installed("irlba")
    k <- 3
    ref <- dense_gplssvd_ref(X0, Y0, MX = NULL, MY = NULL, WX = NULL, WY = NULL,
                             k = k, center = TRUE, scale = TRUE)
    op  <- run_op(X0, Y0, XLW = NULL, YLW = NULL, XRW = NULL, YRW = NULL,
                  k = k, center = TRUE, scale = TRUE, backend = backend)

    # Align signs
    op$u  <- align_signs(ref$u,  op$u)
    op$v  <- align_signs(ref$v,  op$v)
    op$p  <- align_signs(ref$p,  op$p)
    op$q  <- align_signs(ref$q,  op$q)
    op$fi <- align_signs(ref$fi, op$fi)
    op$fj <- align_signs(ref$fj, op$fj)
    op$lx <- align_signs(ref$lx, op$lx)
    op$ly <- align_signs(ref$ly, op$ly)

    testthat::expect_equal(op$d, ref$d, tolerance = 1e-6)
    expect_mats_equal(op$u, ref$u, tol = 1e-6)
    expect_mats_equal(op$v, ref$v, tol = 1e-6)
    expect_mats_equal(op$p,  ref$p,  tol = 1e-6)
    expect_mats_equal(op$q,  ref$q,  tol = 1e-6)
    expect_mats_equal(op$fi, ref$fi, tol = 1e-6)
    expect_mats_equal(op$fj, ref$fj, tol = 1e-6)
    expect_mats_equal(op$lx, ref$lx, tol = 1e-6)
    expect_mats_equal(op$ly, ref$ly, tol = 1e-6)

    d_from_L <- diag(crossprod(op$lx, op$ly))
    testthat::expect_equal(as.numeric(d_from_L), op$d, tolerance = 1e-6)
  })
}

# ---- 2) Diagonal row weights on both X and Y --------------------------------

for (backend in c("eigencore", "irlba")) {
  testthat::test_that(paste("Diagonal row weights match full SVD (backend =", backend, ")"), {
    if (backend == "irlba") testthat::skip_if_not_installed("irlba")
    k <- 3
    set.seed(if (backend == "eigencore") 201 else 202)
    w_row_X <- runif(N)
    w_row_X <- w_row_X / sum(w_row_X)
    w_row_Y <- runif(N)
    w_row_Y <- w_row_Y / sum(w_row_Y)

    ref <- dense_gplssvd_ref(X0, Y0, MX = w_row_X, MY = w_row_Y, WX = NULL, WY = NULL,
                             k = k, center = TRUE, scale = FALSE)
    op  <- run_op(X0, Y0, XLW = Diagonal(x = w_row_X), YLW = Diagonal(x = w_row_Y),
                  XRW = NULL, YRW = NULL, k = k, center = TRUE, scale = FALSE, backend = backend)

    op$u  <- align_signs(ref$u,  op$u)
    op$v  <- align_signs(ref$v,  op$v)
    op$p  <- align_signs(ref$p,  op$p)
    op$q  <- align_signs(ref$q,  op$q)
    op$fi <- align_signs(ref$fi, op$fi)
    op$fj <- align_signs(ref$fj, op$fj)
    op$lx <- align_signs(ref$lx, op$lx)
    op$ly <- align_signs(ref$ly, op$ly)

    testthat::expect_equal(op$d, ref$d, tolerance = 1e-6)
    expect_mats_equal(op$u,  ref$u)
    expect_mats_equal(op$v,  ref$v)
    expect_mats_equal(op$p,  ref$p)
    expect_mats_equal(op$q,  ref$q)
    expect_mats_equal(op$fi, ref$fi)
    expect_mats_equal(op$fj, ref$fj)
    expect_mats_equal(op$lx, ref$lx)
    expect_mats_equal(op$ly, ref$ly)

    d_from_L <- diag(crossprod(op$lx, op$ly))
    testthat::expect_equal(as.numeric(d_from_L), op$d, tolerance = 1e-6)
  })
}

# ---- 3) Full SPD metrics (sparse and dense) on rows and columns --------------

for (backend in c("eigencore", "irlba")) {
  testthat::test_that(paste("Full SPD metrics match full SVD (backend =", backend, ")"), {
    if (backend == "irlba") testthat::skip_if_not_installed("irlba")
    k <- 3

    set.seed(if (backend == "eigencore") 301 else 302)
    MX <- make_spd(N, dense = FALSE) # sparse SPD
    MY <- make_spd(N, dense = TRUE)  # dense SPD
    WX <- make_spd(I, dense = TRUE)
    WY <- make_spd(J, dense = FALSE) # sparse SPD

    ref <- dense_gplssvd_ref(X0, Y0, MX = MX, MY = MY, WX = WX, WY = WY,
                             k = k, center = FALSE, scale = FALSE)
    op  <- run_op(X0, Y0, XLW = MX, YLW = MY, XRW = WX, YRW = WY,
                  k = k, center = FALSE, scale = FALSE, backend = backend)

    op$u  <- align_signs(ref$u,  op$u)
    op$v  <- align_signs(ref$v,  op$v)
    op$p  <- align_signs(ref$p,  op$p)
    op$q  <- align_signs(ref$q,  op$q)
    op$fi <- align_signs(ref$fi, op$fi)
    op$fj <- align_signs(ref$fj, op$fj)
    op$lx <- align_signs(ref$lx, op$lx)
    op$ly <- align_signs(ref$ly, op$ly)

    testthat::expect_equal(op$d, ref$d, tolerance = 1e-6)
    expect_mats_equal(op$u,  ref$u,  tol = 1e-6)
    expect_mats_equal(op$v,  ref$v,  tol = 1e-6)
    expect_mats_equal(op$p,  ref$p,  tol = 1e-6)
    expect_mats_equal(op$q,  ref$q,  tol = 1e-6)
    expect_mats_equal(op$fi, ref$fi, tol = 1e-6)
    expect_mats_equal(op$fj, ref$fj, tol = 1e-6)
    expect_mats_equal(op$lx, ref$lx, tol = 1e-6)
    expect_mats_equal(op$ly, ref$ly, tol = 1e-6)

    d_from_L <- diag(crossprod(op$lx, op$ly))
    testthat::expect_equal(as.numeric(d_from_L), op$d, tolerance = 1e-6)
  })
}

# ---- 4) Sanity: k smaller than rank, and with scaling only on X --------------

for (backend in c("eigencore", "irlba")) {
  testthat::test_that(paste("Partial rank & mixed preprocessing (backend =", backend, ")"), {
    if (backend == "irlba") testthat::skip_if_not_installed("irlba")
    k <- 2
    ref <- dense_gplssvd_ref(X0, Y0, MX = NULL, MY = NULL, WX = NULL, WY = NULL,
                             k = k, center = TRUE, scale = FALSE)
    op  <- run_op(X0, Y0, XLW = NULL, YLW = NULL, XRW = NULL, YRW = NULL,
                  k = k, center = TRUE, scale = FALSE, backend = backend)

    op$u  <- align_signs(ref$u,  op$u)
    op$v  <- align_signs(ref$v,  op$v)
    op$p  <- align_signs(ref$p,  op$p)
    op$q  <- align_signs(ref$q,  op$q)
    op$fi <- align_signs(ref$fi, op$fi)
    op$fj <- align_signs(ref$fj, op$fj)
    op$lx <- align_signs(ref$lx, op$lx)
    op$ly <- align_signs(ref$ly, op$ly)

    testthat::expect_equal(op$d, ref$d, tolerance = 1e-6)
    expect_mats_equal(op$u,  ref$u,  tol = 1e-6)
    expect_mats_equal(op$v,  ref$v,  tol = 1e-6)
    expect_mats_equal(op$p,  ref$p,  tol = 1e-6)
    expect_mats_equal(op$q,  ref$q,  tol = 1e-6)
    expect_mats_equal(op$fi, ref$fi, tol = 1e-6)
    expect_mats_equal(op$fj, ref$fj, tol = 1e-6)
    expect_mats_equal(op$lx, ref$lx, tol = 1e-6)
    expect_mats_equal(op$ly, ref$ly, tol = 1e-6)

    d_from_L <- diag(crossprod(op$lx, op$ly))
    testthat::expect_equal(as.numeric(d_from_L), op$d, tolerance = 1e-6)
  })
}
