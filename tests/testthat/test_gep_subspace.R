library(testthat)
library(Matrix)

# Check the generalized eigen equation S1 V = S2 V diag(lambda) and the
# S2-orthonormality V' S2 V = I. The solver is iterative, so the residual of
# the eigen equation is checked at 1e-3 (relative; observed ~1e-5) and the
# orthonormality at 1e-6 (observed ~1e-8).
check_gep <- function(S1, S2, V, lambda, tol_eq = 1e-3, tol_orth = 1e-6) {
  left <- as.matrix(S1 %*% V)
  right <- as.matrix(S2 %*% (V %*% diag(lambda)))
  expect_equal(left, right, tolerance = tol_eq)
  gram <- as.matrix(t(V) %*% S2 %*% V)
  expect_equal(gram, diag(length(lambda)), tolerance = tol_orth)
}

# Small SPD pencil with a symmetric-whitened reference: the eigenvalues of
# S2^{-1} S1 are those of L^{-1} S1 L^{-T} where S2 = L L' (all real, computed
# via a symmetric eigendecomposition; no complex/Re() issues).
make_gep_problem <- function() {
  set.seed(123)
  n <- 8
  A <- matrix(rnorm(n * n), n, n)
  S2 <- crossprod(A) + diag(n) * 0.1
  B <- matrix(rnorm(n * n), n, n)
  S1 <- crossprod(B) + diag(n) * 0.1
  L <- t(chol(S2))
  Sw <- solve(L, t(solve(L, S1)))
  ev <- eigen((Sw + t(Sw)) / 2, symmetric = TRUE, only.values = TRUE)$values
  list(S1 = Matrix(S1), S2 = Matrix(S2), values = ev, q = 3)
}

test_that("solve_gep_subspace matches direct solution for largest eigenvalues", {
  pb <- make_gep_problem()
  res <- genpca:::solve_gep_subspace(pb$S1, pb$S2, q = pb$q, which = "largest",
                                     reg_S = 1e-6, reg_T = 1e-8)
  expect_equal(res$values, pb$values[seq_len(pb$q)], tolerance = 1e-5)
  check_gep(pb$S1, pb$S2, res$vectors, res$values)
})

test_that("solve_gep_subspace matches direct solution for smallest eigenvalues", {
  pb <- make_gep_problem()
  res <- genpca:::solve_gep_subspace(pb$S1, pb$S2, q = pb$q, which = "smallest",
                                     reg_S = 1e-6, reg_T = 1e-8)
  expect_equal(res$values, rev(pb$values)[seq_len(pb$q)], tolerance = 1e-5)
  check_gep(pb$S1, pb$S2, res$vectors, res$values)
})
