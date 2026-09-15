library(testthat)
library(Matrix)

make_sfpca <- function(seed = 7, n = 60, p = 30, K = 2) {
  set.seed(seed)
  u1 <- sin(seq(0, 2 * pi, length.out = n))
  v1 <- c(rnorm(8), rep(0, p - 8))
  X <- 8 * tcrossprod(u1 / sqrt(sum(u1^2)), v1 / sqrt(sum(v1^2))) +
    matrix(rnorm(n * p, sd = 0.3), n, p)
  spat <- matrix(runif(2 * p), nrow = 2)
  list(fit = sfpca(X, K = K, spat_cds = spat), X = X, n = n, p = p, K = K)
}

test_that("sfpca returns a multivarious bi_projector", {
  s <- make_sfpca()
  fit <- s$fit
  expect_s3_class(fit, "sfpca")
  expect_s3_class(fit, "bi_projector")
  expect_s3_class(fit, "projector")

  expect_equal(multivarious::ncomp(fit), s$K)
  expect_equal(dim(multivarious::scores(fit)), c(s$n, s$K))
  expect_equal(dim(multivarious::components(fit)), c(s$p, s$K))
  expect_length(multivarious::sdev(fit), s$K)
})

test_that("scores equal U D and reconstruct gives the rank-K approximation", {
  s <- make_sfpca()
  fit <- s$fit
  U <- fit$ou
  D <- multivarious::sdev(fit)
  V <- multivarious::components(fit)

  # scores = U D
  expect_equal(multivarious::scores(fit), sweep(U, 2, D, `*`),
               ignore_attr = TRUE)

  # reconstruct = U D V' = scores V'
  rec <- multivarious::reconstruct(fit)
  expect_equal(dim(rec), dim(s$X))
  expect_equal(as.matrix(rec),
               multivarious::scores(fit) %*% t(V), ignore_attr = TRUE)
})

test_that("reconstruct uses t(V) (U D V'), NOT pinv(V), for non-orthogonal loadings", {
  # Two overlapping sparse blocks => materially non-orthogonal loadings.
  # This is the case where scores %*% t(V) and scores %*% pinv(V) differ, and
  # where the inherited reconstruct.bi_projector (pseudoinverse) is WRONG.
  set.seed(99)
  n <- 80
  p <- 40
  ua <- rnorm(n); ub <- rnorm(n)
  va <- c(rnorm(12), rep(0, p - 12))            # support 1:12
  vb <- c(rep(0, 6), rnorm(12), rep(0, p - 18)) # support 7:18 (overlaps va)
  X <- 6 * tcrossprod(ua / sqrt(sum(ua^2)), va / sqrt(sum(va^2))) +
       4 * tcrossprod(ub / sqrt(sum(ub^2)), vb / sqrt(sum(vb^2))) +
       matrix(rnorm(n * p, sd = 0.2), n, p)
  spat <- matrix(runif(2 * p), nrow = 2)

  fit <- sfpca(X, K = 2, spat_cds = spat,
               lambda_u = 0.05, lambda_v = 0.05,
               alpha_u = 0.3, alpha_v = 0.3)

  V <- multivarious::components(fit)
  S <- multivarious::scores(fit)
  # Guard the test's own premise: both components are alive and non-orthogonal,
  # so pinv(V) and t(V) genuinely differ (otherwise the test proves nothing).
  expect_equal(multivarious::ncomp(fit), 2)
  expect_true(all(multivarious::sdev(fit) > 1e-6))
  VtV <- crossprod(V)
  expect_gt(max(abs(VtV - diag(diag(VtV)))), 0.05)

  rec <- as.matrix(multivarious::reconstruct(fit))
  udv  <- S %*% t(V)                          # the sfpca deflation model U D V'
  # pinv(V) for full-column-rank V = (V'V)^{-1} V'; the inherited (wrong) path
  udpi <- S %*% (solve(crossprod(V)) %*% t(V))
  expect_equal(rec, udv, ignore_attr = TRUE, tolerance = 1e-10)
  expect_gt(norm(rec - udpi, "F"), 0.1)       # and it is NOT the pinv version
})

test_that("reconstruct honors comp/rowind/colind subsetting", {
  s <- make_sfpca(seed = 3, K = 2)
  fit <- s$fit
  S <- multivarious::scores(fit)
  V <- multivarious::components(fit)

  r1 <- as.matrix(multivarious::reconstruct(fit, comp = 1))
  expect_equal(r1, S[, 1, drop = FALSE] %*% t(V[, 1, drop = FALSE]),
               ignore_attr = TRUE)

  r_sub <- as.matrix(multivarious::reconstruct(fit, rowind = 1:5, colind = 1:4))
  expect_equal(dim(r_sub), c(5, 4))
  expect_equal(r_sub, S[1:5, , drop = FALSE] %*% t(V[1:4, , drop = FALSE]),
               ignore_attr = TRUE)

  r_null <- as.matrix(multivarious::reconstruct(fit, rowind = NULL, colind = NULL))
  expect_equal(dim(r_null), dim(s$X))
  expect_equal(r_null, S %*% t(V), ignore_attr = TRUE)
})

test_that("penalty parameters are stored as first-class fields", {
  s <- make_sfpca()
  fit <- s$fit
  for (nm in c("lambda_u", "lambda_v", "alpha_u", "alpha_v")) {
    expect_length(fit[[nm]], s$K)
    expect_true(all(is.finite(fit[[nm]])))
  }
  # Reading a stored field via $ must not warn
  expect_silent(fit$lambda_v)
  expect_silent(fit$ov)
})

test_that("deprecated $d and $u still read but warn; $v is the native field", {
  s <- make_sfpca()
  fit <- s$fit

  expect_warning(d_old <- fit$d, "deprecated")
  expect_warning(u_old <- fit$u, "deprecated")
  expect_equal(d_old, multivarious::sdev(fit))
  expect_equal(u_old, fit$ou)

  # $v is the canonical loadings field, identical to components(), no warning
  expect_silent(v_old <- fit$v)
  expect_equal(v_old, multivarious::components(fit))
})

test_that("print.sfpca summarizes the fit without error", {
  s <- make_sfpca()
  expect_output(print(s$fit), "Sparse Functional PCA")
  expect_output(print(s$fit), "components: 2")
})

test_that("K = 1 returns well-formed matrices (no diag(scalar) trap)", {
  s <- make_sfpca(K = 1)
  fit <- s$fit
  expect_equal(dim(multivarious::scores(fit)), c(s$n, 1))
  expect_equal(dim(multivarious::components(fit)), c(s$p, 1))
  expect_length(multivarious::sdev(fit), 1)
})
