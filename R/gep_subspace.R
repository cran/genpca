

# Factor a matrix once with regularization
#' @keywords internal
#' @noRd
factor_mat <- function(M, reg = 1e-3, max_tries = 5) {
 d <- nrow(M)
  for (i in seq_len(max_tries)) {
    M_reg <- M + Matrix::Diagonal(d, reg)
    ch <- try(Matrix::Cholesky(M_reg, LDL = FALSE), silent = TRUE)
    if (!inherits(ch, "try-error")) {
      return(list(ch = ch, final_reg = reg))
    }
    reg <- reg * 10
  }
  stop("Unable to factor matrix even after multiple attempts.")
}


#' @keywords internal
#' @noRd
orthonormalize <- function(X, metric = NULL, reg = 1e-8, max_tries = 12) {
  k <- ncol(X)
  if (k == 0L) {
    return(X[, 0, drop = FALSE])
  }

  if (is.null(metric)) {
    # thin QR; avoids allocating the full Q when d >> q
    QR <- qr(X)
    return(qr.Q(QR, complete = FALSE))
  }

  # Metric orthonormalization: X_hat = X * R^{-1}, where
  # X^T metric X = R^T R.
  G <- crossprod(X, metric %*% X)
  G <- (G + t(G)) / 2

  reg_now <- reg
  R <- NULL
  for (i in seq_len(max_tries)) {
    R <- try(chol(as.matrix(G) + diag(reg_now, k)), silent = TRUE)
    if (!inherits(R, "try-error")) break
    reg_now <- if (reg_now > 0) reg_now * 10 else 1e-10
  }
  if (inherits(R, "try-error")) {
    # Last-resort spectral shift for numerically indefinite small Gram matrices.
    eval <- eigen(as.matrix(G), symmetric = TRUE, only.values = TRUE)$values
    shift <- max(reg_now, -min(eval) + max(reg, 1e-10))
    for (j in seq_len(6)) {
      R <- try(chol(as.matrix(G) + diag(shift, k)), silent = TRUE)
      if (!inherits(R, "try-error")) break
      shift <- shift * 10
    }
  }
  if (inherits(R, "try-error")) {
    stop("Unable to metric-orthonormalize basis: chol failed after regularization and spectral shift.")
  }

  X %*% backsolve(R, diag(1, k))
}


#' Fast sub-space solver for a small block of generalized eigen-pairs
#'
#' Uses pre-conditioned sub-space iteration on the operator
#' \eqn{S_2^{-1} S_1} (or its inverse) to obtain the `q`
#' largest or smallest generalized eigen-values/vectors of
#' \eqn{S_1 v = \lambda S_2 v}.
#'
#' @param S1,S2 Symmetric positive-(semi)definite `dgCMatrix`
#'   (or dense) matrices of the same dimension \eqn{d\times d}.
#' @param q Number of eigen-pairs required (`q << d`).
#' @param which `"largest"` or `"smallest"`.
#' @param max_iter,tol Stopping rule - iteration stops when
#'   `max(abs(lambda_new - lambda_old)/abs(lambda_old)) < tol`.
#' @param V0 Optional `d x q` initial block (will be orthonormalised).
#' @param seed Optional integer seed for reproducible random initialisation.
#' @param reg_S,reg_T Ridge terms added to `S1`/`S2` and the small
#'   `q x q` Gram matrix to guarantee invertibility.
#' @param verbose Logical - print convergence info.
#'
#' @return A list with components
#'   \describe{
#'     \item{values}{length-`q` numeric vector of Ritz eigen-values.}
#'     \item{vectors}{`d x q` matrix, columns are orthonormal eigen-vectors
#'                    in the *original* S-inner-product.}
#'   }
#' @keywords internal
solve_gep_subspace <- function(S1, S2, q = 2, which = c("largest", "smallest"),
                               max_iter = 100, tol = 1e-6, V0 = NULL, seed = NULL,
                               reg_S = 1e-3, reg_T = 1e-6, verbose = FALSE) {
  which <- match.arg(which)
  d <- nrow(S1)

  # Choose operator (avoid explicit inverse)
  # If largest: we iterate with S2^-1 S1 (i.e. solve S2 * V_hat = S1 * V)
  # If smallest: we iterate with S1^-1 S2 (i.e. solve S1 * V_hat = S2 * V)

  if (which == "largest") {
    # Factor S2 once
    s_fact <- factor_mat(S2, reg = reg_S)
    # S2^{-1} S1 V via triangular solves on the Cholesky of S2
    solve_step <- function(V) solve(s_fact$ch, S1 %*% V)
    final_reg_S <- s_fact$final_reg # Store the final reg used
  } else {
    # smallest
    s_fact <- factor_mat(S1, reg = reg_S)
    # S1^{-1} S2 V via triangular solves on the Cholesky of S1
    solve_step <- function(V) solve(s_fact$ch, S2 %*% V)
    final_reg_S <- s_fact$final_reg # Store the final reg used
  }

  if (is.null(V0)) {
    if (!is.null(seed)) set.seed(seed)
    V <- matrix(rnorm(d * q), d, q)
  } else {
    V <- V0
  }

  # Keep iterates S2-orthonormal to make the projected metric matrix stable.
  V <- orthonormalize(V, metric = S2, reg = reg_T)

  lambda_old <- NULL

  for (iter in seq_len(max_iter)) {
    Y <- solve_step(V)

    # Orthonormalize Y in the S2 metric
    Y <- orthonormalize(Y, metric = S2, reg = reg_T)

    # Form S and T efficiently (avoid d x d intermediate products)
    S_mat <- t(Y) %*% (S1 %*% Y)  # q x q
    T_mat <- t(Y) %*% (S2 %*% Y)  # q x q (always S2 here)

    # Symmetrize numerically
    S_mat <- Matrix::forceSymmetric(S_mat)
    T_mat <- Matrix::forceSymmetric(T_mat)

    # Robust Cholesky-whitening of T (small q x q matrix)
    reg <- reg_T
    R <- NULL
    qq <- ncol(T_mat)
    for (tries in 0:7) {
      T_reg <- as.matrix(T_mat) + diag(reg, qq)
      R <- try(chol(T_reg), silent = TRUE)
      if (!inherits(R, "try-error")) break
      reg <- reg * 10
    }
    if (inherits(R, "try-error")) {
      eval_T <- eigen(as.matrix(T_mat), symmetric = TRUE, only.values = TRUE)$values
      shift <- max(reg, -min(eval_T) + max(reg_T, 1e-10))
      for (j in seq_len(6)) {
        R <- try(chol(as.matrix(T_mat) + diag(shift, qq)), silent = TRUE)
        if (!inherits(R, "try-error")) break
        shift <- shift * 10
      }
    }
    if (inherits(R, "try-error")) stop("Unable to chol() the T matrix even after regularization and spectral shift.")

    # C = R^{-T} S R^{-1} (symmetric)
    C <- backsolve(R, t(backsolve(R, t(as.matrix(S_mat)), transpose = TRUE)))
    C <- (C + t(C)) / 2  # Symmetrize

    # Eigen decomposition of whitened matrix
    E <- eigen(C, symmetric = TRUE)

    # Select eigenvalues/vectors based on which
    ord <- if (which == "largest") {
      order(E$values, decreasing = TRUE)
    } else {
      order(E$values, decreasing = FALSE)
    }
    ord <- ord[seq_len(q)]

    lambda <- E$values[ord]
    U <- E$vectors[, ord, drop = FALSE]

    # Transform back: W = R^{-1} U
    W <- backsolve(R, U)

    V_new <- Y %*% W

    if (!is.null(lambda_old)) {
      rel_change <- max(abs(lambda - lambda_old) / pmax(abs(lambda), 1e-12))
      if (rel_change < tol) {
        V <- V_new
        if (verbose) message("Converged in ", iter, " iterations.")
        break
      }
    }

    V <- V_new
    lambda_old <- lambda
  }

  if (iter == max_iter) {
    warning("Reached max_iter without convergence (delta_lambda > tol)")
  }

  # Enforce S2-orthonormality again (with robust regularization).
  V <- orthonormalize(V, metric = S2, reg = reg_T)

  # Recompute lambda as Rayleigh quotients in the now S2-orthonormal basis
  lambda <- diag(crossprod(V, S1 %*% V))

  list(values = as.numeric(lambda), vectors = V, final_reg_S = final_reg_S)
}
