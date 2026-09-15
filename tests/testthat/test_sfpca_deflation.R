library(testthat)
library(Matrix)

test_that("Implicit deflation operators match explicit residual products", {
  set.seed(31)
  n <- 40
  p <- 25
  k <- 3
  X <- Matrix::Matrix(matrix(rnorm(n * p), n, p))
  U <- qr.Q(qr(matrix(rnorm(n * k), n, k)))
  V <- qr.Q(qr(matrix(rnorm(p * k), p, k)))
  d <- c(5, 3, 1)

  X_res <- as.matrix(X) - U %*% (d * t(V))
  ops <- genpca:::sfpca_make_ops(X, U, d, V)

  expect_equal(ops$n, n)
  expect_equal(ops$p, p)
  for (i in 1:5) {
    v <- rnorm(p)
    u <- rnorm(n)
    expect_equal(ops$mv(v), as.numeric(X_res %*% v), tolerance = 1e-12)
    expect_equal(ops$tmv(u), as.numeric(crossprod(X_res, u)), tolerance = 1e-12)
  }

  # Without deflation the operators are plain X products
  ops0 <- genpca:::sfpca_make_ops(X)
  v <- rnorm(p)
  expect_equal(ops0$mv(v), as.numeric(X %*% v), tolerance = 1e-12)
})

test_that("Implicit deflation reproduces explicit-residual extraction", {
  set.seed(77)
  n <- 80
  p <- 50
  K <- 3
  X <- matrix(rnorm(n * p), n, p)
  X[1:40, 1:25] <- X[1:40, 1:25] + 2  # low-rank-ish structure
  spat <- matrix(runif(2 * p), nrow = 2)

  lambda <- 0.1
  alpha <- 0.5
  res_impl <- sfpca(X, K = K, spat_cds = spat,
                    lambda_u = lambda, lambda_v = lambda,
                    alpha_u = alpha, alpha_v = alpha,
                    max_iter = 200, tol = 1e-9)
  impl_u <- res_impl$ou
  impl_v <- multivarious::components(res_impl)
  impl_d <- multivarious::sdev(res_impl)

  # Reference: same per-component problem solved on the explicitly deflated
  # matrix, with initializers computed the same way as sfpca() computes them.
  Omega_u <- Matrix::crossprod(genpca:::second_diff_matrix(n))
  Omega_v <- genpca:::construct_spatial_penalty(spat, k = min(6, p - 1))
  XM <- Matrix::Matrix(X, sparse = FALSE)
  X_res <- XM
  for (k in seq_len(K)) {
    U_prev <- if (k > 1) impl_u[, seq_len(k - 1), drop = FALSE] else NULL
    V_prev <- if (k > 1) impl_v[, seq_len(k - 1), drop = FALSE] else NULL
    d_prev <- if (k > 1) impl_d[seq_len(k - 1)] else numeric(0)
    sv <- genpca:::svd1_deflated(XM, U_prev, d_prev, V_prev)
    ref_k <- genpca:::sfpca_rank1(X_res,
                                  lambda_u = lambda, lambda_v = lambda,
                                  alpha_u = alpha, alpha_v = alpha,
                                  Omega_u = Omega_u, Omega_v = Omega_v,
                                  penalty_u = "l1", penalty_v = "l1",
                                  max_iter = 200, tol = 1e-9, verbose = FALSE,
                                  u_init = as.numeric(sv$u),
                                  v_init = as.numeric(sv$v),
                                  d_init = sv$d[1])
    expect_equal(ref_k$d, impl_d[k], tolerance = 1e-6)
    # u and v are determined only up to a joint sign flip; align before comparing
    sgn <- sign(sum(ref_k$u * impl_u[, k]))
    if (sgn == 0) sgn <- 1
    expect_lt(max(abs(sgn * ref_k$u - impl_u[, k])), 1e-5)
    expect_lt(max(abs(sgn * ref_k$v - impl_v[, k])), 1e-5)
    X_res <- X_res - ref_k$d * Matrix::tcrossprod(ref_k$u, ref_k$v)
  }
})

test_that("Sparse input stays sparse and matches the dense computation", {
  set.seed(13)
  n <- 200
  p <- 100
  X_sp <- Matrix::rsparsematrix(n, p, density = 0.05)
  # add a planted sparse rank-1 signal
  u1 <- rep(0, n); u1[1:50] <- rnorm(50)
  v1 <- rep(0, p); v1[1:20] <- rnorm(20)
  X_sp <- X_sp + 3 * Matrix::Matrix(tcrossprod(u1, v1), sparse = TRUE)
  X_sp <- Matrix::drop0(X_sp)
  spat <- matrix(runif(2 * p), nrow = 2)

  res_sp <- sfpca(X_sp, K = 2, spat_cds = spat,
                  lambda_u = 0.05, lambda_v = 0.05,
                  alpha_u = 0.5, alpha_v = 0.5, tol = 1e-9)
  res_dn <- sfpca(as.matrix(X_sp), K = 2, spat_cds = spat,
                  lambda_u = 0.05, lambda_v = 0.05,
                  alpha_u = 0.5, alpha_v = 0.5, tol = 1e-9)

  sp_u <- res_sp$ou; dn_u <- res_dn$ou
  sp_v <- multivarious::components(res_sp); dn_v <- multivarious::components(res_dn)
  expect_equal(multivarious::sdev(res_sp), multivarious::sdev(res_dn),
               tolerance = 1e-6)
  for (k in 1:2) {
    # align signs before comparing
    s <- sign(sum(sp_u[, k] * dn_u[, k]))
    expect_lt(max(abs(s * sp_u[, k] - dn_u[, k])), 1e-5)
    expect_lt(max(abs(s * sp_v[, k] - dn_v[, k])), 1e-5)
  }
})
