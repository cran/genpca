library(testthat)
library(Matrix)

test_that("lambda_max = ||b||_inf is the exact zero-solution threshold", {
  set.seed(11)
  p <- 40
  b <- rnorm(p)
  Omega <- Matrix::crossprod(genpca:::second_diff_matrix(p))
  S <- genpca:::as_dgc(Matrix::Diagonal(p) + 0.5 * Omega)
  lam_max <- max(abs(b))

  sol_at <- genpca:::sfpca_cd_solve(S, b, numeric(p), lam_max, "l1")
  expect_true(all(sol_at$x == 0))

  sol_below <- genpca:::sfpca_cd_solve(S, b, numeric(p), 0.95 * lam_max, "l1")
  expect_gt(sum(sol_below$x != 0), 0)
})

test_that("sfpca_res_fnorm2 matches the explicit residual Frobenius norm", {
  set.seed(21)
  n <- 30
  p <- 20
  k <- 3
  X <- matrix(rnorm(n * p), n, p)
  U <- matrix(rnorm(n * k), n, k)
  V <- matrix(rnorm(p * k), p, k)
  d <- c(4, 2, 1)

  F2 <- genpca:::sfpca_res_fnorm2(sum(X^2), X, U, d, V)
  F2_explicit <- sum((X - U %*% (d * t(V)))^2)
  expect_equal(F2, F2_explicit, tolerance = 1e-10)

  # No deflation: passes through
  expect_equal(genpca:::sfpca_res_fnorm2(sum(X^2), X), sum(X^2))
})

test_that("default_alpha is scale-free and bounds the subproblem conditioning", {
  set.seed(31)
  p <- 50
  spat <- matrix(runif(2 * p), nrow = 2)
  Omega <- genpca:::construct_spatial_penalty(spat, k = 5)

  a1 <- genpca:::default_alpha(Omega)
  a2 <- genpca:::default_alpha(10 * Omega)
  # alpha * Omega is invariant to rescaling of Omega
  expect_equal(a1, 10 * a2, tolerance = 1e-8)

  # cond(I + alpha * Omega) <= 2: largest eigenvalue of alpha*Omega is 1
  ev <- eigen(as.matrix(Matrix::Diagonal(p) + a1 * Omega),
              symmetric = TRUE, only.values = TRUE)$values
  expect_lte(max(ev), 2 + 1e-8)
  expect_gte(min(ev), 1 - 1e-8)

  # Zero penalty gives alpha = 0, not Inf
  expect_equal(genpca:::default_alpha(Matrix::Matrix(0, p, p, sparse = TRUE)), 0)
})

test_that("BIC path selection prefers a sparse fit on planted sparse data", {
  set.seed(41)
  n <- 120
  p <- 80
  u1 <- rnorm(n)
  v1 <- rep(0, p)
  v1[1:10] <- rnorm(10, sd = 2)
  X <- tcrossprod(u1, v1) + matrix(rnorm(n * p, sd = 0.5), n, p)

  ops <- genpca:::sfpca_make_ops(X)
  sv <- genpca:::svd1_deflated(X)
  S_v <- genpca:::as_dgc(Matrix::Diagonal(p))
  u_fix <- as.numeric(sv$u)
  b_v <- ops$tmv(u_fix)

  sel <- genpca:::sfpca_select_lambda(b_v, S_v, sum(X^2), n * p,
                                      sqrt(sum(u_fix^2)), "l1",
                                      nlambda = 10, lambda_min_ratio = 1e-2)
  expect_length(sel$lambdas, 10)
  expect_length(sel$bic, 10)
  expect_equal(sel$lambda, sel$lambdas[sel$index])
  # First path point is lambda_max (null model), which must not win here
  expect_equal(sel$lambdas[1], max(abs(b_v)))
  expect_gt(sel$index, 1)

  # The selected lambda produces a solution supported (mostly) on the truth
  sol <- genpca:::sfpca_cd_solve(S_v, b_v, numeric(p), sel$lambda, "l1")
  supp <- which(sol$x != 0)
  expect_lt(length(supp), p / 2)
  expect_gte(length(intersect(supp, 1:10)), 8)
})

test_that("sfpca defaults (BIC lambda, scale-free alpha) recover planted structure", {
  set.seed(51)
  n <- 100
  p <- 60
  spat <- matrix(runif(2 * p), nrow = 2)
  u1 <- sin(seq(0, 2 * pi, length.out = n))
  v1 <- rep(0, p)
  v1[1:12] <- rnorm(12, sd = 2)
  # Signal strong enough that the BIC RSS gain justifies a dense u
  # (log(np)/np per df sets the bar; a marginal signal correctly loses)
  X <- 10 * tcrossprod(u1 / sqrt(sum(u1^2)), v1 / sqrt(sum(v1^2))) +
    matrix(rnorm(n * p, sd = 0.2), n, p)

  res <- sfpca(X, K = 1, spat_cds = spat)

  expect_gt(multivarious::sdev(res)[1], 0)
  # Selected penalties are recorded and positive
  expect_gt(res$lambda_u[1], 0)
  expect_gt(res$lambda_v[1], 0)
  expect_gt(res$alpha_u[1], 0)
  expect_gt(res$alpha_v[1], 0)

  # Component aligns with the planted signal
  expect_gt(abs(cor(res$ou[, 1], u1)), 0.8)
  expect_gt(abs(cor(multivarious::components(res)[, 1], v1)), 0.8)
})

test_that("uthresh/vthresh are deprecated with a warning and ignored", {
  set.seed(61)
  n <- 40
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  spat <- matrix(runif(2 * p), nrow = 2)

  expect_warning(
    res_dep <- sfpca(X, K = 1, spat_cds = spat, uthresh = 0.9),
    "deprecated"
  )
  res_def <- sfpca(X, K = 1, spat_cds = spat)
  expect_equal(multivarious::sdev(res_dep), multivarious::sdev(res_def))
  expect_equal(res_dep$lambda_v, res_def$lambda_v)
})

test_that("nlambda and lambda_min_ratio are validated", {
  set.seed(134)
  X <- matrix(rnorm(200), 20, 10)
  spat <- matrix(runif(20), nrow = 2)
  expect_error(sfpca(X, K = 1, spat_cds = spat, nlambda = 1), "nlambda")
  expect_error(sfpca(X, K = 1, spat_cds = spat, lambda_min_ratio = 1.5),
               "lambda_min_ratio")
})

test_that("negative or non-finite alpha is rejected (keeps S = I + alpha Omega SPD)", {
  set.seed(142)
  X <- matrix(rnorm(200), 20, 10)
  spat <- matrix(runif(20), nrow = 2)
  expect_error(sfpca(X, K = 1, spat_cds = spat, alpha_v = -1), "alpha_v")
  expect_error(sfpca(X, K = 1, spat_cds = spat, alpha_u = Inf), "alpha_u")
})
