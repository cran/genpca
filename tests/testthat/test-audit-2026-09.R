# Regression tests for the 2026-09 audit (plans/2026-09-04-audit-remediation.md).
# Each block names the audit item it covers. Blocks whose phase has not landed
# yet are skipped; remove the skip when the phase lands.

audit_X <- function(n = 60, p = 8, seed = 1) {
  set.seed(seed)
  matrix(rnorm(n * p), n, p)
}
audit_spd <- function(p, seed = 2) {
  set.seed(seed)
  crossprod(matrix(rnorm(p * p), p)) / p + diag(p) * 0.1
}

# ---- 8.1: metric-first truncation ------------------------------------------
test_that("8.1: maxeig never truncates the metric (exact or refuse)", {
  set.seed(1)
  n <- 50; p <- 6; eps <- 1e-2
  Qo <- qr.Q(qr(matrix(rnorm(p * p), p)))
  lam <- c(1, 0.9, 0.8, 0.7, 0.6, eps)
  A <- Qo %*% diag(lam) %*% t(Qo); A <- (A + t(A)) / 2
  Z <- scale(matrix(rnorm(n * p), n), scale = FALSE)
  X <- Z %*% diag(c(1, 1, 1, 1, 1, 1 / eps)) %*% t(Qo)
  Ah <- Qo %*% diag(sqrt(lam)) %*% t(Qo)
  true_sdev <- sqrt(eigen(Ah %*% crossprod(X) %*% Ah, symmetric = TRUE)$values[1:3])

  f_full <- genpca(X, A = A, ncomp = 3, method = "eigen", preproc = multivarious::pass(), maxeig = Inf)
  expect_equal(f_full$sdev, true_sdev, tolerance = 1e-6)
  # A positive definite metric larger than maxeig is factored exactly.
  f_small <- genpca(X, A = A, ncomp = 3, method = "eigen", preproc = multivarious::pass(), maxeig = 4)
  expect_equal(f_small$sdev, true_sdev, tolerance = 1e-6)
  # A singular general metric larger than maxeig is refused, never approximated.
  A_sing <- Qo %*% diag(c(1, 0.9, 0.8, 0.7, 0.6, 0)) %*% t(Qo); A_sing <- (A_sing + t(A_sing)) / 2
  expect_error(
    genpca(X, A = A_sing, ncomp = 3, method = "eigen", preproc = multivarious::pass(), maxeig = 4),
    "maxeig"
  )
  f_sing <- genpca(X, A = A_sing, ncomp = 3, method = "eigen", preproc = multivarious::pass(), maxeig = Inf)
  expect_equal(length(f_sing$sdev), 3L)
})

# ---- 8.2: PSD vs PD semantics, tolerance forwarding, weight vectors --------
test_that("8.2: is_psd/is_pd are distinct and ensure_spd returns PD", {
  Z <- matrix(0, 3, 3)
  expect_true(is_psd(Z))
  expect_false(is_pd(Z))
  r <- ensure_spd(Z)
  expect_true(min(eigen(as.matrix(r), symmetric = TRUE)$values) > 0)
  M <- diag(c(1, -5e-7))
  expect_false(is_psd(M))
  expect_true(min(eigen(as.matrix(ensure_spd(M, tol = 1e-12)), symmetric = TRUE)$values) > 0)
})

test_that("8.2: negative diagonal weights are treated identically by every backend", {
  X <- audit_X(40, 5)
  w_bad <- c(1, 1, 1, 1, -5e-7)
  w_tiny <- c(1, 1, 1, 1, -1e-12)
  ref <- NULL
  for (m in c("eigen", "spectra", "deflation", "randomized")) {
    expect_error(genpca(X, A = w_bad, ncomp = 2, method = m, constraints_remedy = "error"),
                 "non-negative", label = m)
    f <- suppressMessages(genpca(X, A = w_tiny, ncomp = 2, method = m, seed_randomized = 1L))
    if (is.null(ref)) ref <- f$sdev
    expect_equal(f$sdev, ref, tolerance = 1e-6, label = m)
  }
})

# ---- 8.3: symmetry is measured, not assumed; clip really clips -------------
test_that("8.3: asymmetric metrics error under every remedy; tiny asymmetry is averaged", {
  X <- audit_X(40, 2)
  A_asym <- matrix(c(2, 0.9, 0.1, 2), 2)          # relative asymmetry ~0.28
  for (rem in c("error", "ridge", "clip", "identity")) {
    expect_error(genpca(X, A = A_asym, ncomp = 1, constraints_remedy = rem), "symmetric", label = rem)
  }
  A_sym <- matrix(c(2, 0.5, 0.5, 2), 2)
  A_near <- A_sym; A_near[1, 2] <- A_near[1, 2] + 1e-13
  f_near <- genpca(X, A = A_near, ncomp = 1)
  f_sym <- genpca(X, A = A_sym, ncomp = 1)
  expect_equal(f_near$sdev, f_sym$sdev, tolerance = 1e-10)
  C <- crossprod(X); Ca <- C; Ca[1, 2] <- C[1, 2] + 1
  expect_error(genpca_cov(Ca, ncomp = 1), "symmetric")
})

test_that("8.3: clip_psd output has no negative eigenvalues", {
  r <- clip_psd(diag(c(1, -5e-7)))
  expect_true(min(eigen(as.matrix(r), symmetric = TRUE)$values) >= 0)
})

# ---- 8.4: scale equivariance of rank decisions -----------------------------
test_that("8.4: every backend is invariant to rescaling X (primal and dual, down and up)", {
  shapes <- list(primal = c(60, 8), dual = c(8, 60))
  for (sh in names(shapes)) {
    n <- shapes[[sh]][1]; p <- shapes[[sh]][2]
    X <- audit_X(n, p)
    Am <- audit_spd(p); Mm <- diag(runif(n, 0.5, 1.5))
    for (m in c("eigen", "spectra", "deflation", "randomized")) {
      ref <- genpca(X, A = Am, M = Mm, ncomp = 4, method = m, preproc = multivarious::pass(), seed_randomized = 1L)
      for (cc in c(1e-7, 1e-5, 1e-3, 1e3, 1e6, 1e9)) {
        f <- genpca(X * cc, A = Am, M = Mm, ncomp = 4, method = m, preproc = multivarious::pass(), seed_randomized = 1L)
        lab <- paste(sh, m, cc)
        expect_equal(length(f$sdev), length(ref$sdev), label = lab)
        expect_equal(f$sdev / cc, ref$sdev, tolerance = 1e-6, label = lab)
        # loadings must not be zeroed at any scale (dual-branch guard regression)
        expect_true(all(colSums(abs(as.matrix(multivarious::components(f)))) > 0), label = lab)
      }
    }
  }
})

# ---- 8.7: gpca_mle reports the objective at the returned metrics -----------
test_that("8.7: gpca_mle loglik equals the penalized objective at the returned (M, A, fit)", {
  X <- audit_X(30, 6, seed = 2)
  n <- nrow(X); p <- ncol(X); lambda <- 1e-3
  pen_obj <- function(res) {
    E <- X - multivarious::reconstruct(res$fit)
    Sr <- Matrix::solve(res$M); Sc <- Matrix::solve(res$A)
    quad <- sum(E * as.matrix(res$M %*% E %*% res$A))
    pen <- p * lambda * sum(Matrix::diag(res$M)) + n * lambda * sum(Matrix::diag(res$A))
    ldr <- as.numeric(Matrix::determinant(Sr, logarithm = TRUE)$modulus)
    ldc <- as.numeric(Matrix::determinant(Sc, logarithm = TRUE)$modulus)
    -0.5 * (p * ldr + n * ldc + quad + pen)
  }
  for (sf in c("none", "trace", "det")) {
    res <- suppressWarnings(gpca_mle(X, ncomp = 2, max_iter = 6, lambda = lambda, scale_fix = sf))
    expect_equal(res$loglik, pen_obj(res), tolerance = 1e-6, label = sf)
  }
})

# ---- 8.8: metric repair is explicit ----------------------------------------
test_that("8.8: an indefinite metric errors by default and warns when repaired", {
  X <- audit_X(40, 2)
  A_ind <- matrix(c(1, 2, 2, 1), 2)
  expect_error(genpca(X, A = A_ind, ncomp = 1), "positive semi-definite")
  expect_warning(
    f <- genpca(X, A = A_ind, ncomp = 1, constraints_remedy = "ridge"),
    class = "genpca_metric_repaired"
  )
  expect_true(min(eigen(as.matrix(f$A), symmetric = TRUE)$values) > 0)
  rep <- repair_metric(A_ind, method = "ridge")
  expect_s3_class(attr(rep, "repair_report"), "metric_repair_report")
})
