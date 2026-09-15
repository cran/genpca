# Phase 3 hardening: pin recently fixed behaviors and core invariants
# (scores/project semantics, relative deflation thresholds, randomized RNG
# hygiene, clip vs ridge remediation, cache keying, colind reconstruction).

library(testthat)
library(Matrix)

test_that("gpca_mle penalized loglik path is monotone and metrics are SPD", {
  # Every block update (GPCA rank-k fit, Sigma_r, Sigma_c) exactly minimizes
  # the shared penalized objective, so loglik_path must be non-decreasing up
  # to numerical noise.
  set.seed(11)
  X <- matrix(rnorm(30 * 12), 30, 12)
  for (sf in c("trace", "det", "none")) {
    res <- suppressWarnings(
      gpca_mle(X, ncomp = 3, max_iter = 8, scale_fix = sf, verbose = FALSE))
    path <- res$loglik_path
    expect_true(all(is.finite(path)))
    expect_lte(length(path), 8)
    # Account separately for reciprocal rescaling and final reevaluation.
    if (sf == "none") {
      expect_equal(res$loglik, tail(path, 1), tolerance = 1e-8)
      expect_equal(res$loglik_rescale_delta, 0, tolerance = 1e-8)
    } else {
      expect_equal(res$loglik, tail(path, 1) + res$loglik_rescale_delta + res$loglik_refit_delta)
    }
    if (length(path) > 1) {
      slack <- 1e-8 * pmax(abs(path[-length(path)]), 1)
      expect_true(all(diff(path) >= -slack),
                  info = paste("non-monotone path for scale_fix =", sf))
    }
    expect_equal(dim(res$A), c(12, 12))
    expect_equal(dim(res$M), c(30, 30))
    expect_true(genpca:::is_spd(res$A))
    expect_true(genpca:::is_spd(res$M))
    # returned fit was computed with the returned (canonicalized) metrics
    expect_equal(as.matrix(res$fit$M), as.matrix(res$M), tolerance = 1e-12)
    expect_equal(as.matrix(res$fit$A), as.matrix(res$A), tolerance = 1e-12)
  }

  # Low-rank data: previously the worst oscillation case
  set.seed(12)
  L <- matrix(rnorm(30 * 2), 30, 2) %*% t(matrix(rnorm(12 * 2), 12, 2))
  X2 <- L + 0.05 * matrix(rnorm(30 * 12), 30, 12)
  res2 <- suppressWarnings(
    gpca_mle(X2, ncomp = 2, max_iter = 10, scale_fix = "none", verbose = FALSE))
  path2 <- res2$loglik_path
  if (length(path2) > 1) {
    slack2 <- 1e-8 * pmax(abs(path2[-length(path2)]), 1)
    expect_true(all(diff(path2) >= -slack2))
  }
})

test_that("rpls with lambda=0 matches the PLS-SVD reference and predicts training Y", {
  set.seed(61)
  n <- 60; p <- 12; q <- 6; K <- 3
  Tn <- matrix(rnorm(n * K), n, K)
  X <- Tn %*% matrix(rnorm(K * p), K, p) + 0.01 * matrix(rnorm(n * p), n, p)
  Y <- Tn %*% matrix(rnorm(K * q), K, q) + 0.01 * matrix(rnorm(n * q), n, q)

  fit <- rpls(X, Y, K = K, lambda = 0, penalty = "ridge",
              preproc_x = multivarious::center(),
              preproc_y = multivarious::center())
  Xc <- scale(X, scale = FALSE)
  Yc <- scale(Y, scale = FALSE)
  sv <- svd(crossprod(Xc, Yc))

  # Unpenalized first component equals the leading PLS-SVD pair
  expect_gt(abs(cor(fit$vx[, 1], sv$u[, 1])), 1 - 1e-6)
  expect_gt(abs(cor(fit$vy[, 1], sv$v[, 1])), 1 - 1e-6)
  expect_gt(abs(cor(Xc %*% fit$vx[, 1], Xc %*% sv$u[, 1])), 1 - 1e-6)

  # Later components span (essentially) the same latent subspace:
  # cosines of the principal angles between the two K-dim latent spaces
  cosines <- svd(crossprod(qr.Q(qr(Xc %*% fit$vx)),
                           qr.Q(qr(Xc %*% sv$u[, 1:K]))))$d
  expect_true(all(cosines > 0.99))

  # Predictive sanity: transferred X correlates strongly with training Y
  ty <- multivarious::transfer(fit, X, from = "X", to = "Y")
  expect_gt(cor(as.vector(as.matrix(ty)), as.vector(Y)), 0.85)
})

test_that("truncate.genpca slices every component-indexed field consistently", {
  set.seed(63)
  n <- 25; p <- 10; k <- 2
  X <- matrix(rnorm(n * p), n, p)
  Mm <- crossprod(matrix(rnorm(n * n), n, n)) + diag(n) * 0.5
  Aa <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p) * 0.5
  fit <- genpca(X, A = Aa, M = Mm, ncomp = 5, preproc = multivarious::center())
  tr <- multivarious::truncate(fit, k)

  expect_s3_class(tr, "genpca")
  expect_equal(multivarious::ncomp(tr), k)
  expect_equal(as.matrix(tr$v), as.matrix(fit$v[, 1:k]), ignore_attr = TRUE)
  expect_equal(as.matrix(tr$s), as.matrix(fit$s[, 1:k]), ignore_attr = TRUE)
  expect_equal(multivarious::sdev(tr), multivarious::sdev(fit)[1:k])
  expect_equal(as.matrix(tr$ov), as.matrix(fit$ov[, 1:k]), ignore_attr = TRUE)
  expect_equal(as.matrix(tr$ou), as.matrix(fit$ou[, 1:k]), ignore_attr = TRUE)
  expect_equal(as.matrix(tr$u), as.matrix(fit$u[, 1:k]), ignore_attr = TRUE)
  expect_equal(tr$propv, fit$propv[1:k])
  expect_equal(tr$cumv, fit$cumv[1:k])

  expect_equal(reconstruct(tr), reconstruct(fit, comp = 1:k), tolerance = 1e-12)
})

test_that("transfer.cross_projector: legacy source/target args match from/to", {
  set.seed(64)
  n <- 40; p <- 8; q <- 5; K <- 2
  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1:3] %*% matrix(rnorm(3 * q), 3, q) + 0.01 * matrix(rnorm(n * q), n, q)
  fit <- rpls(X, Y, K = K, lambda = 0, penalty = "ridge")

  t_new <- multivarious::transfer(fit, X, from = "X", to = "Y")
  t_old <- multivarious::transfer(fit, X, source = "X", target = "Y")
  expect_identical(t_new, t_old)

  t_new_yx <- multivarious::transfer(fit, Y, from = "Y", to = "X")
  t_old_yx <- multivarious::transfer(fit, Y, source = "Y", target = "X")
  expect_identical(t_new_yx, t_old_yx)
})

test_that("transfer.cross_projector round trip recovers a noiseless rank-K system", {
  set.seed(65)
  n <- 60; p <- 12; q <- 6; K <- 3
  Tn <- matrix(rnorm(n * K), n, K)
  X <- Tn %*% matrix(rnorm(K * p), K, p)          # exactly rank K
  Y <- X %*% matrix(rnorm(p * q), p, q)           # noiseless linear map
  fit <- rpls(X, Y, K = K, lambda = 0, penalty = "ridge")

  xy <- multivarious::transfer(fit, X, from = "X", to = "Y")
  back <- multivarious::transfer(fit, xy, from = "Y", to = "X")
  # Round trip goes through two least-squares solves; accumulated error is
  # BLAS-dependent (observed 2e-6 on win-builder R-devel vs 1e-9 locally).
  expect_lt(norm(as.matrix(back) - X, "F") / norm(X, "F"), 1e-5)

  # X transferred to Y-space correlates strongly with the true Y
  expect_gt(cor(as.vector(as.matrix(xy)), as.vector(Y)), 0.85)
})

test_that("method='randomized' is deterministic given seed_randomized and leaves .Random.seed alone", {
  set.seed(41)
  X <- matrix(rnorm(60 * 30), 60, 30)

  set.seed(1)
  saved <- get(".Random.seed", envir = globalenv())
  fitA <- genpca(X, ncomp = 4, method = "randomized", seed_randomized = 7,
                 preproc = multivarious::pass())
  expect_identical(get(".Random.seed", envir = globalenv()), saved)

  fitB <- genpca(X, ncomp = 4, method = "randomized", seed_randomized = 7,
                 preproc = multivarious::pass())
  fitC <- genpca(X, ncomp = 4, method = "randomized", seed_randomized = 8,
                 preproc = multivarious::pass())

  expect_identical(fitA$sdev, fitB$sdev)
  expect_identical(as.matrix(fitA$ov), as.matrix(fitB$ov))
  expect_false(isTRUE(all.equal(as.matrix(fitA$ov), as.matrix(fitC$ov))))
})

test_that("scores/project coherence under the ou-D convention with dense SPD metrics", {
  set.seed(51)
  n <- 25; p <- 10; k <- 5
  X <- matrix(rnorm(n * p), n, p)
  Mm <- crossprod(matrix(rnorm(n * n), n, n)) + diag(n) * 0.5
  Aa <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p) * 0.5
  fit <- genpca(X, A = Aa, M = Mm, ncomp = k, preproc = multivarious::center())

  sc <- as.matrix(multivarious::scores(fit))
  pr <- as.matrix(multivarious::project(fit, X))
  expect_equal(sc, pr, tolerance = 1e-8, ignore_attr = TRUE)

  # scores = ou D (the paper's generalized PCs)
  expect_equal(sc, sweep(as.matrix(fit$ou), 2, multivarious::sdev(fit), `*`),
               tolerance = 1e-8, ignore_attr = TRUE)

  # components v = A ov; u = M ou
  expect_equal(as.matrix(fit$v), as.matrix(Aa %*% fit$ov),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(as.matrix(fit$u), as.matrix(Mm %*% fit$ou),
               tolerance = 1e-8, ignore_attr = TRUE)

  # metric orthonormality
  expect_equal(as.matrix(crossprod(fit$ou, Mm %*% fit$ou)), diag(k),
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(as.matrix(crossprod(fit$ov, Aa %*% fit$ov)), diag(k),
               tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("reconstruct with colind honors centering (column subset matches X)", {
  set.seed(52)
  n <- 25; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  fit <- genpca(X, ncomp = p, preproc = multivarious::center())
  rc <- reconstruct(fit, colind = c(2, 4))
  expect_equal(as.matrix(rc), X[, c(2, 4)], tolerance = 1e-8,
               ignore_attr = TRUE)
})

test_that("constraints_remedy 'clip' differs from 'ridge' and preserves the positive spectrum", {
  set.seed(21)
  n <- 12; p <- 8
  X <- matrix(rnorm(n * p), n, p)
  G <- crossprod(matrix(rnorm(p * p), p, p))
  ee <- eigen(G, symmetric = TRUE)
  vals <- ee$values
  vals[p] <- -0.5
  vals[p - 1] <- -0.1
  A_ind <- ee$vectors %*% (vals * t(ee$vectors))
  A_ind <- (A_ind + t(A_ind)) / 2  # indefinite metric

  fit_clip <- genpca(X, A = A_ind, ncomp = 3, constraints_remedy = "clip",
                     preproc = multivarious::pass())
  fit_ridge <- genpca(X, A = A_ind, ncomp = 3, constraints_remedy = "ridge",
                      preproc = multivarious::pass())

  # A real spectral clip is not a diagonal ridge shift
  expect_false(isTRUE(all.equal(fit_clip$sdev, fit_ridge$sdev)))

  eig_clip <- eigen(as.matrix(fit_clip$A), symmetric = TRUE,
                    only.values = TRUE)$values
  expect_gte(min(eig_clip), -1e-8)
  # largest eigenvalues of the clipped metric equal the original positive ones
  pos <- sort(vals[vals > 0], decreasing = TRUE)
  expect_equal(eig_clip[1:3], pos[1:3], tolerance = 1e-8)
})

test_that("deflation thresholds are scale invariant", {
  set.seed(31)
  X <- matrix(rnorm(40 * 15), 40, 15)

  set.seed(99)
  f1 <- genpca(X, ncomp = 5, method = "deflation",
               preproc = multivarious::center())
  set.seed(99)
  f2 <- genpca(X * 1e-8, ncomp = 5, method = "deflation",
               preproc = multivarious::center())

  expect_equal(multivarious::ncomp(f1), 5)
  expect_equal(multivarious::ncomp(f2), multivarious::ncomp(f1))
  expect_equal(f2$sdev / 1e-8, f1$sdev, tolerance = 1e-4)
})

test_that("sfpca_cd_solve_cpp rejects non-finite inputs", {
  S <- Matrix::sparseMatrix(i = 1:3, j = 1:3, x = c(1, 2, 3))
  ok <- genpca:::sfpca_cd_solve_cpp(S, c(1, 0, 2), rep(0, 3), 0.1, 0L, 3.7,
                                    100L, 1e-8)
  expect_true(all(is.finite(ok$x)))
  expect_error(
    genpca:::sfpca_cd_solve_cpp(S, c(1, NA, 2), rep(0, 3), 0.1, 0L, 3.7,
                                100L, 1e-8),
    "finite")
  expect_error(
    genpca:::sfpca_cd_solve_cpp(S, c(1, Inf, 2), rep(0, 3), 0.1, 0L, 3.7,
                                100L, 1e-8),
    "finite")
  expect_error(
    genpca:::sfpca_cd_solve_cpp(S, c(1, 0, 2), c(0, NaN, 0), 0.1, 0L, 3.7,
                                100L, 1e-8),
    "finite")
})

test_that("gmd cache keys distinguish nearly identical small-magnitude matrices", {
  genpca::gmd_clear_cache()
  on.exit(genpca::gmd_clear_cache(), add = TRUE)

  A1 <- diag(1e-9, 3)
  A2 <- diag(2e-9, 3)
  L1 <- genpca:::get_chol_lower_dense(A1)
  L2 <- genpca:::get_chol_lower_dense(A2)

  expect_equal(L1 %*% t(L1), A1, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(L2 %*% t(L2), A2, tolerance = 1e-12, ignore_attr = TRUE)
  expect_false(isTRUE(all.equal(L1, L2)))
})
