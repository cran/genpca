library(testthat)
library(Matrix)

# Helpers -------------------------------------------------------------------

# Sparse SPD S = I + alpha * (second-difference penalty)
make_S <- function(n, alpha) {
  D <- genpca:::second_diff_matrix(n)
  genpca:::as_dgc(Matrix::Diagonal(n) + alpha * Matrix::crossprod(D))
}

subproblem_objective <- function(x, S, b, lambda, penalty) {
  0.5 * as.numeric(Matrix::crossprod(x, S %*% x)) - sum(b * x) +
    genpca:::penalty_value(x, penalty, lambda)
}

# Independent R-side KKT residual for the l1 subproblem
kkt_residual_l1 <- function(x, S, b, lambda) {
  g <- as.numeric(S %*% x - b)
  r_nz <- if (any(x != 0)) max(abs(g[x != 0] + lambda * sign(x[x != 0]))) else 0
  r_z <- if (any(x == 0)) max(0, max(abs(g[x == 0])) - lambda) else 0
  max(r_nz, r_z)
}

# Reference solver: long proximal-gradient (ISTA) run
ista_reference <- function(S, b, lambda, n_iter = 5000) {
  L <- max(abs(eigencore::eigs_sym(S, 1, which = "LM")$values))
  x <- rep(0, length(b))
  for (i in seq_len(n_iter)) {
    y <- x - as.numeric(S %*% x - b) / L
    x <- sign(y) * pmax(0, abs(y) - lambda / L)
  }
  x
}

# Tests ----------------------------------------------------------------------

test_that("CD solver satisfies KKT conditions for the l1 subproblem", {
  set.seed(42)
  n <- 80
  S <- make_S(n, alpha = 2.5)
  b <- rnorm(n, sd = 3)

  for (lambda in c(0, 0.1, 1, 5)) {
    sol <- genpca:::sfpca_cd_solve(S, b, rep(0, n), lambda, "l1")
    expect_lt(kkt_residual_l1(sol$x, S, b, lambda), 1e-6)
  }

  # lambda = 0 reduces to the linear system S x = b
  sol0 <- genpca:::sfpca_cd_solve(S, b, rep(0, n), 0, "l1")
  x_direct <- as.numeric(Matrix::solve(S, b))
  expect_lt(max(abs(sol0$x - x_direct)), 1e-6)

  # lambda >= ||b||_inf gives the zero solution
  sol_big <- genpca:::sfpca_cd_solve(S, b, rep(0, n), max(abs(b)) + 1, "l1")
  expect_true(all(sol_big$x == 0))
})

test_that("CD solution matches a long proximal-gradient reference (l1)", {
  set.seed(7)
  n <- 60
  S <- make_S(n, alpha = 1.0)
  b <- rnorm(n, sd = 2)
  lambda <- 0.5

  x_cd <- genpca:::sfpca_cd_solve(S, b, rep(0, n), lambda, "l1")$x
  x_ref <- ista_reference(S, b, lambda)

  f_cd <- subproblem_objective(x_cd, S, b, lambda, "l1")
  f_ref <- subproblem_objective(x_ref, S, b, lambda, "l1")
  expect_lte(f_cd, f_ref + 1e-8)
  expect_lt(max(abs(x_cd - x_ref)), 1e-4)
})

test_that("Warm starts reproduce the cold-start solution", {
  set.seed(11)
  n <- 50
  S <- make_S(n, alpha = 0.5)
  b <- rnorm(n)
  lambda <- 0.3

  cold <- genpca:::sfpca_cd_solve(S, b, rep(0, n), lambda, "l1")$x
  warm <- genpca:::sfpca_cd_solve(S, b, cold + rnorm(n, sd = 0.1), lambda, "l1")$x
  expect_lt(max(abs(cold - warm)), 1e-6)
})

test_that("SCAD coordinate updates are exact univariate minimizers", {
  # With diagonal S the subproblem separates; compare each coordinate against
  # a brute-force univariate minimization.
  set.seed(3)
  n <- 30
  s_diag <- runif(n, 1, 3)
  S <- genpca:::as_dgc(Matrix::Diagonal(x = s_diag))
  b <- rnorm(n, sd = 4)
  lambda <- 1
  a <- 3.7

  sol <- genpca:::sfpca_cd_solve(S, b, rep(0, n), lambda, "scad")$x

  scad_pen <- function(x) genpca:::penalty_value(x, "scad", lambda, a)
  for (j in seq_len(n)) {
    h <- function(x) 0.5 * s_diag[j] * x^2 - b[j] * x + scad_pen(x)
    # Brute-force minimizer over a fine grid (covers all SCAD regimes)
    grid <- seq(-2 * abs(b[j]) / s_diag[j] - 1, 2 * abs(b[j]) / s_diag[j] + 1,
                length.out = 20001)
    x_bf <- grid[which.min(vapply(grid, h, numeric(1)))]
    expect_lt(h(sol[j]) - h(x_bf), 1e-6)
  }
})

test_that("Rank-1 objective trace is monotone non-decreasing (l1)", {
  set.seed(21)
  n <- 60
  p <- 40
  # Smooth-in-time u, smooth-in-space v, plus noise
  u_true <- sin(seq(0, 2 * pi, length.out = n))
  spat <- matrix(seq(0, 1, length.out = p), nrow = 1)
  v_true <- cos(2 * pi * spat[1, ])
  X <- 5 * tcrossprod(u_true, v_true) + matrix(rnorm(n * p, sd = 0.5), n, p)

  Omega_u <- Matrix::crossprod(genpca:::second_diff_matrix(n))
  Omega_v <- genpca:::construct_spatial_penalty(spat, k = 4)

  res <- genpca:::sfpca_rank1(Matrix::Matrix(X),
                              lambda_u = 0.2, lambda_v = 0.2,
                              alpha_u = 1, alpha_v = 1,
                              Omega_u = Omega_u, Omega_v = Omega_v,
                              penalty_u = "l1", penalty_v = "l1",
                              max_iter = 100, tol = 1e-9, verbose = FALSE,
                              exact_inner = TRUE)

  expect_gt(length(res$obj_trace), 1)
  expect_true(all(diff(res$obj_trace) >= -1e-7))
  # Subproblems solved to optimality at the last iteration
  expect_lt(res$kkt_u, 1e-5)
  expect_lt(res$kkt_v, 1e-5)
  # And the smooth structure is recovered
  expect_gt(abs(cor(res$u, u_true)), 0.95)
  expect_gt(abs(cor(res$v, v_true)), 0.95)
})

test_that("Overwhelming penalty yields a clean zero component", {
  set.seed(5)
  n <- 30
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  Omega_u <- Matrix::crossprod(genpca:::second_diff_matrix(n))
  spat <- matrix(runif(2 * p), nrow = 2)
  Omega_v <- genpca:::construct_spatial_penalty(spat, k = 4)

  res <- genpca:::sfpca_rank1(Matrix::Matrix(X),
                              lambda_u = 1e6, lambda_v = 1e6,
                              alpha_u = 0.1, alpha_v = 0.1,
                              Omega_u = Omega_u, Omega_v = Omega_v,
                              penalty_u = "l1", penalty_v = "l1",
                              max_iter = 50, tol = 1e-8, verbose = FALSE)
  expect_equal(res$d, 0)
  expect_true(all(res$u == 0))
  expect_true(all(res$v == 0))
})
