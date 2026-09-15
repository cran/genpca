#' Internal Utilities for Partial Eigen, Adaptive Rank, and Sqrt Transforms
#'
#' These functions are used by both \code{genpls} and \code{genplscorr} to handle
#' partial eigen expansions, diagonal/identity shortcuts, adaptive rank selection,
#' and row/column transformations for data embeddings.
#'
#' @name plsutils
#' @keywords internal
#' @importFrom Matrix diag isDiagonal crossprod
NULL

#' Check if constraint matrix is identity or purely diagonal with all != 1
#'
#' @param M A matrix (often `dsCMatrix`) or `NULL`.
#' @param eps Numeric tolerance
#' @return TRUE if \code{M} is a (Matrix-based) diagonal with all or partial diag
#' @keywords internal
is_identity_or_diag <- function(M, eps = 1e-15) {
  (inherits(M, "Matrix") && Matrix::isDiagonal(M))
}

#' Partial eigen decomposition up to rank k
#'
#' @param M symmetric PSD matrix
#' @param k maximum rank
#' @param which "LA" for largest algebraic
#' @param eps small numeric for safety
#' @param tol tolerance for the iterative eigensolver
#' @return list(Q=..., lam=...) with length(lam)=k_eff
#' @keywords internal
#' @noRd
partial_eig_once <- function(M, k = 50, which = "LA", eps = 1e-15, tol = 1e-6) {
  # clamp k
  k_eff <- min(k, nrow(M), ncol(M))
  k_eff <- max(k_eff, 1)

  # For small or full-rank requests, base eigen() is more accurate and avoids
  # asking the iterative solver for k >= n eigenpairs.
  if (nrow(M) <= 50 || k_eff >= nrow(M)) {
    es <- eigen(as.matrix(M), symmetric = TRUE)
    lam <- pmax(es$values[seq_len(k_eff)], 0)
    Q   <- es$vectors[, seq_len(k_eff), drop = FALSE]
  } else {
    es <- .top_eigs_sym(M, k_eff, which, tol = tol)
    lam <- pmax(es$values, 0)
    Q   <- es$vectors
  }

  if (!is.matrix(Q) || nrow(Q) != nrow(M)) {
    stop("partial_eig_once: dimension mismatch or no eigenvectors returned.")
  }

  list(Q = Q, lam = lam)
}

#' Decide adaptive rank to capture var_threshold fraction of total variance
#'
#' @param M PSD matrix
#' @param var_threshold fraction (0..1) of variance to capture
#' @param max_k maximum rank
#' @inheritParams partial_eig_once
#' @return list(Q, lam) truncated to minimal r s.t cumsum(lam)/sum(lam)>=var_threshold
#' @keywords internal
#' @noRd
decide_adaptive_rank <- function(M, which = "LA", eps = 1e-15, tol = 1e-6,
                                 var_threshold = 0.99, max_k = 200)
{
  out <- partial_eig_once(M, k = max_k, which = which, eps = eps, tol = tol)
  lam <- out$lam
  Q   <- out$Q
  cumsums <- cumsum(lam)
  total <- sum(lam)
  r <- which(cumsums >= var_threshold * total)[1]
  if (is.na(r)) {
    r <- length(lam)
  }
  list(
    Q   = Q[, 1:r, drop = FALSE],
    lam = lam[1:r]
  )
}

#' partial_eig_approx: unify user-specified rank or adaptive
#'
#' @param M PSD matrix
#' @param user_rank numeric>0 or NA/NULL/0 => adaptive
#' @param var_threshold fraction
#' @param max_k maximum rank
#' @keywords internal
#' @noRd
partial_eig_approx <- function(M, user_rank,
                               var_threshold = 0.99, max_k = 200,
                               which = "LA", eps = 1e-15, tol = 1e-6)
{
  if (!is.null(user_rank) && !is.na(user_rank) && user_rank > 0) {
    # direct partial
    out <- partial_eig_once(M, k = user_rank, which = which, eps = eps, tol = tol)
    list(Q = out$Q, lams = out$lam)
  } else {
    # adaptive
    adapt <- decide_adaptive_rank(M, which = which, eps = eps, tol = tol,
                                  var_threshold = var_threshold, max_k = max_k)
    list(Q = adapt$Q, lams = adapt$lam)
  }
}

#' Build a closure that multiplies data by sqrt(M)
#'
#' @param M PSD matrix or NULL => identity
#' @param user_rank numeric or NULL => adaptive
#' @inheritParams partial_eig_approx
#' @return function(mat) => mat' = sqrt(M)* mat
#' @keywords internal
#' @noRd
build_sqrt_mult <- function(M, user_rank, var_threshold = 0.99, max_k = 200,
                            which = "LA", eps = 1e-15, tol = 1e-6)
{
  if (is.null(M) || is_identity_or_diag(M)) {
    # Delegate identity/diagonal to shared metric ops for consistency
    ops <- .metric_operators(M, if (!is.null(M)) nrow(M) else NULL)
    return(function(mat) ops$mult_sqrt(mat))
  } else {
    out <- partial_eig_approx(M, user_rank, var_threshold, max_k,
                              which = which, eps = eps, tol = tol)
    Q <- out$Q
    lam <- out$lams
    return(function(mat) {
      alpha <- crossprod(Q, mat)
      for (i in seq_along(lam)) {
        alpha[i, ] <- alpha[i, ] * sqrt(pmax(lam[i], eps))
      }
      Q %*% alpha
    })
  }
}

#' Build a closure that multiplies data by M^(-1/2)
#'
#' @keywords internal
#' @noRd
build_invsqrt_mult <- function(M, user_rank, var_threshold = 0.99, max_k = 200,
                               which = "LA", eps = 1e-15, tol = 1e-6)
{
  if (is.null(M) || is_identity_or_diag(M)) {
    ops <- .metric_operators(M, if (!is.null(M)) nrow(M) else NULL)
    return(function(vec) ops$mult_invsqrt(vec))
  } else {
    out <- partial_eig_approx(M, user_rank, var_threshold, max_k,
                              which = which, eps = eps, tol = tol)
    Q   <- out$Q
    lam <- out$lams
    return(function(vec) {
      alpha <- crossprod(Q, vec)
      for (i in seq_along(lam)) {
        if (lam[i] > eps) {
          alpha[i, ] <- alpha[i, ] / sqrt(lam[i])
        } else {
          alpha[i, ] <- 0
        }
      }
      Q %*% alpha
    })
  }
}

#' row_transform, col_transform helpers
#'
#' @keywords internal
#' @noRd
row_transform <- function(mat, ffun) ffun(mat)

#' @noRd
#'@keywords internal
col_transform <- function(mat, ffun) {
  mt <- t(mat)
  mt2 <- ffun(mt)
  t(mt2)
}
