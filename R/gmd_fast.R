#' Check whether a metric matrix is diagonal
#' @param M a matrix
#' @return logical
#' @keywords internal
is_diagonal_metric <- function(M) {
  inherits(M, "ddiMatrix") || Matrix::isDiagonal(M)
}

diag_metric_values <- function(M) {
  as.numeric(Matrix::diag(M))
}

should_use_topk <- function(k, min_dim, topk, auto_topk, topk_ratio, topk_min_dim) {
  if (!isTRUE(topk)) return(FALSE)
  if (!isTRUE(auto_topk)) return(TRUE)
  if (min_dim < topk_min_dim) return(FALSE)
  (k / min_dim) <= topk_ratio
}

# Components kept after an iterative or dense decomposition: singular values
# above rank_rtol times the largest one (relative, so invariant to rescaling).
.keep_components <- function(d, rank_rtol) {
  dmax <- if (length(d)) max(d, 0) else 0
  is.finite(d) & d > 0 & d > rank_rtol * dmax
}

.empty_gmd <- function(n, p) {
  list(u = matrix(0, n, 0), v = matrix(0, p, 0),
       ou = matrix(0, n, 0), ov = matrix(0, p, 0), d = numeric(0))
}

#' Factor a metric as A = F F' for the whitened-operator SVD
#'
#' Diagonal, dense Cholesky, sparse (CHOLMOD, permuted) Cholesky and
#' PSD-singular (eigen) factors share one closure interface:
#' `apply(V) = F V`, `apply_t(U) = F' U`, `solve_t(Z) = F^{+T} Z` (maps the
#' singular vectors of the whitened operator back to metric-orthonormal
#' factors), `mat` (F itself, for the dense fallback) and `ncol` (columns of
#' F: the dimension for Cholesky factors, the numerical rank for eigen
#' factors). A metric that is not positive definite needs a dense
#' eigendecomposition, which is refused above `dense_maxn` rows.
#' @keywords internal
#' @noRd
.metric_factor <- function(A, metric_rtol = .metric_rtol_default(), cache = TRUE,
                           dense_maxn = 5000L, name = "metric", allow_eigen = TRUE) {
  if (!inherits(A, "Matrix")) A <- Matrix::Matrix(A, sparse = FALSE)
  n <- nrow(A)

  if (is_diagonal_metric(A)) {
    d <- .clamp_weights(diag_metric_values(A), metric_rtol, name)
    s <- sqrt(d)
    # Exact diagonal weights have no eigendecomposition roundoff: retain every
    # positive weight in both the forward factor and its inverse.
    inv <- ifelse(d > 0, 1 / s, 0)
    return(list(kind = "diag", ncol = n, mat = Matrix::Diagonal(n, x = s),
                apply = function(V) s * V,
                apply_t = function(U) s * U,
                solve_t = function(Z) inv * as.matrix(Z)))
  }

  A <- symmetrize_or_stop(A, name = name)

  if (methods::is(A, "sparseMatrix")) {
    ch <- suppressWarnings(tryCatch(
      Matrix::Cholesky(A, LDL = FALSE, perm = TRUE, super = TRUE),
      error = function(e) NULL))
    if (!is.null(ch)) {
      # P A P' = L L'  =>  A = F F' with F = P' L
      L <- if (exists("expand1", envir = asNamespace("Matrix"))) {
        Matrix::expand1(ch, "L")
      } else {
        methods::as(ch, "CsparseMatrix")
      }
      perm <- ch@perm + 1L
      # Same pivot test as the dense branch: a numerically singular factor
      # would amplify null-space noise in solve_t, so fall through to eigen.
      piv <- as.numeric(Matrix::diag(L))
      if (min(piv)^2 > metric_rtol * .metric_scale(A)) {
      Fm <- L[order(perm), , drop = FALSE]
      Lt <- Matrix::t(L)
      return(list(kind = "sparse_chol", ncol = n, mat = Fm,
                  apply = function(V) Fm %*% V,
                  apply_t = function(U) Matrix::crossprod(Fm, U),
                  solve_t = function(Z) {
                    # F^{-T} Z = P' L^{-T} Z
                    W <- as.matrix(Matrix::solve(Lt, as.matrix(Z)))
                    out <- W
                    out[perm, ] <- W
                    out
                  }))
      }
    }
    if (!allow_eigen) return(NULL)
    if (n > dense_maxn) {
      stop("Metric ", name, " is not positive definite and has ", n, " rows (> ", dense_maxn,
           "): a dense eigendecomposition would be needed. Supply a positive definite ",
           "metric, use method = 'deflation', or raise maxeig.", call. = FALSE)
    }
    A <- Matrix::Matrix(as.matrix(A), sparse = FALSE)
  }

  # Dense: Cholesky when positive definite (cached), eigen factor otherwise.
  scale <- .metric_scale(A)
  L <- tryCatch(if (isTRUE(cache)) get_chol_lower_dense(A) else t(chol(as.matrix(A))),
                error = function(e) NULL)
  if (!is.null(L) && min(diag(L))^2 > metric_rtol * scale) {
    Lt <- t(L)
    return(list(kind = "chol", ncol = n, mat = L,
                apply = function(V) L %*% V,
                apply_t = function(U) crossprod(L, U),
                solve_t = function(Z) backsolve(Lt, as.matrix(Z))))
  }
  if (!allow_eigen) return(NULL)
  if (n > dense_maxn) {
    stop("Metric ", name, " is not positive definite and has ", n, " rows (> ", dense_maxn,
         "): a dense eigendecomposition would be needed. Supply a positive definite ",
         "metric, use method = 'deflation', or raise maxeig.", call. = FALSE)
  }
  es <- eigen(as.matrix(A), symmetric = TRUE)
  lam <- es$values
  if (min(lam) < -metric_rtol * scale) {
    stop(name, " must be positive semi-definite (min eigenvalue ", signif(min(lam), 3), ")",
         call. = FALSE)
  }
  keep <- which(lam > metric_rtol * max(lam, 0))
  if (!length(keep)) stop(name, " is (numerically) zero.", call. = FALSE)
  Vk <- es$vectors[, keep, drop = FALSE]
  sk <- sqrt(lam[keep])
  Fm <- Vk * rep(sk, each = nrow(Vk))          # V_k Lambda^{1/2}
  list(kind = "eigen", ncol = length(keep), mat = Fm,
       apply = function(V) Fm %*% V,
       apply_t = function(U) crossprod(Fm, U),
       solve_t = function(Z) Vk %*% (as.matrix(Z) / sk))   # V_k Lambda^{-1/2} Z
}

# Diagonal metrics: partial SVD of the scaled data matrix.
gmd_fast_diag <- function(X, q_diag, r_diag, k, tol, maxit, topk, rank_rtol, metric_rtol) {
  n <- nrow(X)
  p <- ncol(X)
  min_dim <- min(n, p)
  k_use <- min(k, min_dim)

  q_diag <- .clamp_weights(q_diag, metric_rtol, "Q")
  r_diag <- .clamp_weights(r_diag, metric_rtol, "R")
  q_sqrt <- sqrt(q_diag)
  r_sqrt <- sqrt(r_diag)
  q_invsqrt <- ifelse(q_diag > 0, 1 / q_sqrt, 0)
  r_invsqrt <- ifelse(r_diag > 0, 1 / r_sqrt, 0)

  # Target: singular values of Q^{1/2} X R^{1/2}.
  Xw <- if (inherits(X, "Matrix")) {
    (q_sqrt * X) %*% Matrix::Diagonal(p, x = r_sqrt)
  } else {
    sweep(q_sqrt * X, 2, r_sqrt, `*`)
  }

  sv <- NULL
  if (isTRUE(topk) && k_use >= 1 && k_use < min_dim) {
    sv <- tryCatch(.top_svd(Xw, k_use, tol = tol), error = function(e) NULL)
    if (!is.null(sv) && !isTRUE(sv$converged)) sv <- NULL
  }
  if (is.null(sv)) {
    sv_full <- base::svd(as.matrix(Xw), nu = k_use, nv = k_use)
    sv <- list(d = sv_full$d[seq_len(k_use)],
               u = sv_full$u[, seq_len(k_use), drop = FALSE],
               v = sv_full$v[, seq_len(k_use), drop = FALSE])
  }

  d <- as.numeric(sv$d)
  keep <- .keep_components(d, rank_rtol)
  if (!any(keep)) return(.empty_gmd(n, p))
  d <- d[keep]
  Uw <- as.matrix(sv$u)[, keep, drop = FALSE]
  Vw <- as.matrix(sv$v)[, keep, drop = FALSE]

  ov <- r_invsqrt * Vw
  ou <- q_invsqrt * Uw
  components <- r_sqrt * Vw                    # R ov
  scores <- q_diag * as.matrix(X %*% components) # Q ou D
  list(u = scores, v = components, ou = ou, ov = ov, d = d)
}

# General metrics: partial SVD of the whitened operator B = F_Q' X F_R,
# where Q = F_Q F_Q' and R = F_R F_R'. B'B = F_R' X'QX F_R has the GMD
# eigenvalues, ov = F_R^{-T} V and ou = F_Q^{-T} U are the metric-orthonormal
# factors. Neither metric is ever squared-rooted explicitly and only
# products with X, X' and the factors are needed per iteration.
gmd_spectra_general <- function(X, Q, R, k, tol, use_topk, cache, rank_rtol, metric_rtol,
                                dense_maxn = 5000L) {
  n <- nrow(X)
  p <- ncol(X)
  primal <- p <= n
  # The small-side metric may need an eigen factor (bounded by dense_maxn);
  # the large-side metric is only used through products unless it is
  # positive definite, in which case its Cholesky factor gives the SVD form.
  FQ <- .metric_factor(Q, metric_rtol, cache, dense_maxn = dense_maxn, name = "Q",
                       allow_eigen = !primal)
  FR <- .metric_factor(R, metric_rtol, cache, dense_maxn = dense_maxn, name = "R",
                       allow_eigen = primal)

  if (!is.null(FQ) && !is.null(FR)) {
    return(gmd_spectra_svd(X, Q, R, FQ, FR, k, tol, use_topk, rank_rtol))
  }
  gmd_spectra_sym(X, Q, R, if (primal) FR else FQ, primal, k, tol, use_topk, rank_rtol)
}

# SVD form: B = F_Q' X F_R with Q = F_Q F_Q', R = F_R F_R'. B'B = F_R' X'QX F_R
# has the GMD eigenvalues, ov = F_R^{-T} V and ou = F_Q^{-T} U are the
# metric-orthonormal factors. Only products with X, X' and the factors are
# needed per iteration.
gmd_spectra_svd <- function(X, Q, R, FQ, FR, k, tol, use_topk, rank_rtol) {
  n <- nrow(X)
  p <- ncol(X)
  nB <- FQ$ncol
  pB <- FR$ncol
  k_use <- min(k, nB, pB)
  if (k_use < 1) return(.empty_gmd(n, p))

  op <- function(V, args = NULL) FQ$apply_t(X %*% FR$apply(V))
  opt <- function(U, args = NULL) FR$apply_t(Matrix::crossprod(X, FQ$apply(U)))

  sv <- NULL
  if (isTRUE(use_topk) && k_use < min(nB, pB)) {
    sv <- tryCatch(.top_svd(op, k_use, tol = tol, adjoint = opt, dim = c(nB, pB)),
                   error = function(e) NULL)
    if (!is.null(sv) && !isTRUE(sv$converged)) sv <- NULL
  }
  if (is.null(sv)) {
    B <- as.matrix(FQ$apply_t(X %*% FR$mat))   # nB x pB, same footprint as X
    sv_full <- base::svd(B, nu = k_use, nv = k_use)
    sv <- list(d = sv_full$d[seq_len(k_use)],
               u = sv_full$u[, seq_len(k_use), drop = FALSE],
               v = sv_full$v[, seq_len(k_use), drop = FALSE])
  }

  d <- as.numeric(sv$d)
  o <- order(d, decreasing = TRUE)
  d <- d[o]
  U <- as.matrix(sv$u)[, o, drop = FALSE]
  V <- as.matrix(sv$v)[, o, drop = FALSE]
  keep <- .keep_components(d, rank_rtol)
  if (!any(keep)) return(.empty_gmd(n, p))
  d <- d[keep]
  U <- U[, keep, drop = FALSE]
  V <- V[, keep, drop = FALSE]

  ov <- as.matrix(FR$solve_t(V))
  ou <- as.matrix(FQ$solve_t(U))
  components <- as.matrix(R %*% ov)                    # R ov
  scores <- as.matrix(Q %*% ou) * rep(d, each = n)     # Q ou D
  list(u = scores, v = components, ou = ou, ov = ov, d = d)
}

# Symmetric small-side form, used when the large-side metric is singular (so
# it must not be factored): T = F' (X'QX) F (primal, p x p) or F' (X R X') F
# (dual, n x n) with F the small-side factor. Same algebra as gmdLA; the
# large-side metric only appears in products.
gmd_spectra_sym <- function(X, Q, R, Fs, primal, k, tol, use_topk, rank_rtol) {
  n <- nrow(X)
  p <- ncol(X)
  G <- if (primal) Matrix::crossprod(X, Q %*% X) else X %*% (R %*% Matrix::t(X))
  T <- as.matrix(Fs$apply_t(G %*% Fs$mat))
  T <- 0.5 * (T + t(T))
  dim_t <- nrow(T)
  k_use <- min(k, dim_t)
  if (k_use < 1) return(.empty_gmd(n, p))

  eig <- NULL
  if (isTRUE(use_topk) && k_use < dim_t - 1L && dim_t > 500L) {
    eig <- tryCatch(.top_eigs_sym(T, k_use, "LA", tol = 1e-10), error = function(e) NULL)
    if (!is.null(eig) && !isTRUE(eig$converged)) eig <- NULL
  }
  if (is.null(eig)) {
    es <- eigen(T, symmetric = TRUE)
    eig <- list(values = es$values[seq_len(k_use)], vectors = es$vectors[, seq_len(k_use), drop = FALSE])
  }
  lam <- as.numeric(eig$values)
  keep <- lam > 0 & lam > rank_rtol^2 * max(lam, 0)
  if (!any(keep)) return(.empty_gmd(n, p))
  lam <- lam[keep]
  Z <- as.matrix(eig$vectors)[, keep, drop = FALSE]
  d <- sqrt(lam)

  if (primal) {
    ov <- as.matrix(Fs$solve_t(Z))                     # p x k, R-orthonormal
    components <- as.matrix(R %*% ov)
    ou <- as.matrix(X %*% components) / rep(d, each = n) # X R ov = ou D
  } else {
    ou <- as.matrix(Fs$solve_t(Z))                     # n x k, Q-orthonormal
    ov <- as.matrix(Matrix::crossprod(X, Q %*% ou)) / rep(d, each = p) # X'Q ou = ov D
    components <- as.matrix(R %*% ov)
  }
  scores <- as.matrix(Q %*% ou) * rep(d, each = n)
  list(u = scores, v = components, ou = ou, ov = ov, d = d)
}

#' Generalized matrix decomposition via partial SVD of the whitened operator
#'
#' Computes the generalized SVD of X with row metric Q and column metric R,
#' equivalent to the eigendecomposition used by \code{\link{genpca}} with
#' \code{method = "eigen"}. The metrics are factored once (`Q = F_Q F_Q'`,
#' `R = F_R F_R'`; diagonal, dense or sparse Cholesky, or an eigen factor for
#' singular metrics) and the top-k singular triplets of the implicit operator
#' `F_Q' X F_R` are computed with \pkg{eigencore}; a dense SVD is used when
#' few components are not requested or the iterative solver does not
#' converge. `gmd_fast_cpp()` is an alias kept for existing callers.
#'
#' @section When is this fast:
#' \itemize{
#'   \item \code{k << min(n, p)}: only the top-k triplets are computed
#'   \item Repeated calls with the same dense Q or R: Cholesky factors are cached
#'   \item Sparse metrics: only sparse factors and products are formed
#' }
#' A positive definite metric on the big side of X costs one Cholesky of that
#' dimension; a singular one is never factored (symmetric small-side form).
#'
#' @param X numeric matrix (n x p)
#' @param Q,R constraints (weights/metrics) for rows/cols. Must be symmetric
#'   positive (semi-)definite. Can be dense matrices, sparse matrices, or
#'   diagonal matrices.
#' @param k number of components to extract (must be >= 1 and <= min(n, p))
#' @param tol convergence tolerance of the iterative solver. Default 1e-9.
#' @param maxit unused (kept for compatibility).
#' @param seed unused (kept for compatibility); results do not depend on the
#'   R random stream.
#' @param topk logical; use the iterative top-k solver when \code{k < min(n, p)}.
#'   Set to FALSE to force a dense SVD of the whitened operator.
#' @param cache logical; cache dense Cholesky factors across calls.
#'   Defaults to TRUE. Use \code{\link{gmd_clear_cache}} to clear.
#' @param auto_topk logical; when TRUE (default), use top-k only when
#'   \code{k/min(n,p)} is small and \code{min(n,p)} is large enough.
#' @param topk_ratio threshold used by \code{auto_topk}. If
#'   \code{k/min(n,p) <= topk_ratio}, top-k is used. Default 0.08.
#' @param topk_min_dim minimum \code{min(n,p)} required before top-k is used
#'   under \code{auto_topk}. Default 200.
#' @param diag_fast logical; if TRUE (default) and both constraints are
#'   diagonal, use a weighted-SVD fast path.
#' @param rank_rtol relative cutoff on singular values: components with
#'   \code{d_j <= rank_rtol * d_1} are dropped. Default 1e-6.
#' @param metric_rtol relative tolerance for metric validation and null-space
#'   detection. Default \code{sqrt(.Machine$double.eps)}.
#' @param dense_maxn a singular general metric on the small side of X needs a
#'   dense eigendecomposition; refuse it above this many rows (the
#'   \code{maxeig} argument of \code{\link{genpca}}). A singular metric on
#'   the large side is never factored: the solver switches to the symmetric
#'   small-side formulation, in which that metric only appears in products.
#'
#' @return A list with components:
#'   \describe{
#'     \item{u}{n x k matrix of metric-weighted scores \code{Q ou D}}
#'     \item{v}{p x k matrix of components \code{R ov}}
#'     \item{ou,ov}{metric-orthonormal factors}
#'     \item{d}{length-k vector of singular values}
#'     \item{k}{number of components returned (may be < requested if rank-deficient)}
#'   }
#'
#' @seealso \code{\link{genpca}} for the high-level interface,
#'   \code{\link{gmd_clear_cache}} to clear the Cholesky cache
#' @keywords internal
#' @importFrom methods as is
gmd_spectra <- function(X, Q, R, k, tol = 1e-9, maxit = 1000L, seed = 1234L,
                        topk = TRUE, cache = TRUE, auto_topk = TRUE,
                        topk_ratio = 0.08, topk_min_dim = 200L,
                        diag_fast = TRUE, rank_rtol = 1e-6,
                        metric_rtol = .metric_rtol_default(), dense_maxn = 5000L) {
  n <- nrow(X)
  p <- ncol(X)
  if (!is.numeric(k) || length(k) != 1 || k < 1) {
    stop("k must be a single positive integer >= 1")
  }
  k <- as.integer(k)
  if (k > min(n, p)) {
    warning("k (", k, ") exceeds min(n, p) = ", min(n, p), "; will return at most ", min(n, p), " components")
    k <- min(n, p)
  }
  if (!is.numeric(maxit) || length(maxit) != 1 || maxit < 1 || maxit != floor(maxit)) {
    stop("maxit must be a single positive integer >= 1")
  }
  if (!is.numeric(topk_ratio) || length(topk_ratio) != 1 || topk_ratio <= 0 || topk_ratio > 1) {
    stop("topk_ratio must be a single number in (0, 1].")
  }
  if (!is.numeric(topk_min_dim) || length(topk_min_dim) != 1 || topk_min_dim < 2 || topk_min_dim != floor(topk_min_dim)) {
    stop("topk_min_dim must be a single integer >= 2.")
  }
  if (!is.numeric(rank_rtol) || length(rank_rtol) != 1 || !is.finite(rank_rtol) || rank_rtol < 0) {
    stop("rank_rtol must be a single non-negative number.")
  }
  topk_min_dim <- as.integer(topk_min_dim)

  if (!inherits(Q, "Matrix")) Q <- Matrix::Matrix(Q, sparse = FALSE)
  if (!inherits(R, "Matrix")) R <- Matrix::Matrix(R, sparse = FALSE)

  use_topk <- should_use_topk(k, min(n, p), topk, auto_topk, topk_ratio, topk_min_dim)

  res <- if (isTRUE(diag_fast) && is_diagonal_metric(Q) && is_diagonal_metric(R)) {
    gmd_fast_diag(X, diag_metric_values(Q), diag_metric_values(R), k,
                  tol = tol, maxit = maxit, topk = use_topk,
                  rank_rtol = rank_rtol, metric_rtol = metric_rtol)
  } else {
    gmd_spectra_general(X, Q, R, k, tol = tol, use_topk = use_topk, cache = cache,
                        rank_rtol = rank_rtol, metric_rtol = metric_rtol,
                        dense_maxn = dense_maxn)
  }
  res$d <- as.vector(res$d)
  res$k <- length(res$d)
  res
}

#' @rdname gmd_spectra
#' @keywords internal
gmd_fast_cpp <- gmd_spectra

# Metric orthonormalization for small block matrices. `jitter` and `tol` are
# relative to the scale of the Gram matrix, so the result is invariant to
# rescaling A.
metric_orthonormalize <- function(A, applyM, jitter = 1e-10, tol = 1e-12) {
  if (ncol(A) == 0) return(A)

  MA <- applyM(A)
  G <- Matrix::crossprod(A, MA)
  G <- 0.5 * (G + t(G))
  G <- as.matrix(G)
  gscale <- max(diag(G), 0)
  C <- tryCatch(chol(G + diag(max(jitter, 0) * gscale, nrow(G))),
                error = function(e) NULL)
  if (!is.null(C)) {
    B <- as.matrix(A %*% solve(C))
    check <- as.matrix(Matrix::crossprod(B, applyM(B)))
    if (all(is.finite(check)) && norm(check - diag(ncol(B)), "I") <= 1e-8)
      return(B %*% solve(chol(check)))
  }

  # A regularized Gram cannot certify rank or metric orthogonality. Build
  # an independent Euclidean basis first, then whiten its small metric Gram.
  scaled <- as.matrix(A)
  scales <- sqrt(colSums(scaled^2))
  scaled <- sweep(scaled, 2, ifelse(scales > 0, scales, 1), "/")
  sv <- svd(scaled, nv = 0)
  eps <- .Machine$double.eps
  independent <- sv$d > eps * max(dim(A)) * max(sv$d, 0) & sv$d > 0
  if (!any(independent)) return(matrix(0, nrow(A), 0))
  basis <- sv$u[, independent, drop = FALSE]
  G <- as.matrix(crossprod(basis, applyM(basis)))
  eg <- eigen(0.5 * (G + t(G)), symmetric = TRUE)
  keep <- eg$values > eps * nrow(G) * max(eg$values, 0) & eg$values > 0
  if (!any(keep)) return(matrix(0, nrow(A), 0))
  basis %*% sweep(eg$vectors[, keep, drop = FALSE], 2,
                  sqrt(eg$values[keep]), "/")
}

random_sign_matrix <- function(nrow, ncol, seed = NULL) {
  if (!is.null(seed)) {
    # Never clobber the caller's RNG stream (CRAN policy): save and restore
    # .Random.seed around the seeded draw.
    old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      get(".Random.seed", envir = globalenv(), inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
          rm(".Random.seed", envir = globalenv())
        }
      } else {
        assign(".Random.seed", old_seed, envir = globalenv())
      }
    }, add = TRUE)
    set.seed(seed)
  }
  matrix(sample(c(-1, 1), size = nrow * ncol, replace = TRUE), nrow = nrow, ncol = ncol)
}

random_normal_matrix <- function(nrow, ncol, seed = NULL) {
  if (!is.null(seed)) {
    # Never clobber the caller's RNG stream (CRAN policy): save and restore
    # .Random.seed around the seeded draw.
    old_seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      get(".Random.seed", envir = globalenv(), inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
          rm(".Random.seed", envir = globalenv())
        }
      } else {
        assign(".Random.seed", old_seed, envir = globalenv())
      }
    }, add = TRUE)
    set.seed(seed)
  }
  matrix(rnorm(nrow * ncol), nrow = nrow, ncol = ncol)
}

gmd_randomized_polish <- function(X, applyQ, applyR, U, V, iters = 1L, jitter = 1e-10,
                                  polish_tol = 0) {
  if (iters <= 0L || ncol(U) == 0) {
    return(list(u = U, v = V, d = numeric(ncol(U))))
  }

  d <- numeric(ncol(U))
  d_prev <- NULL
  for (it in seq_len(iters)) {
    Y <- X %*% applyR(V)
    U <- metric_orthonormalize(Y, applyQ, jitter = jitter)
    Z <- Matrix::crossprod(X, applyQ(U))
    V <- metric_orthonormalize(Z, applyR, jitter = jitter)
    RV <- applyR(V)
    # T = U^T Q X R V; reuse Z = X^T Q U to avoid another pass through X.
    T <- Matrix::crossprod(Z, RV)
    sv <- svd(as.matrix(T))
    U <- U %*% sv$u
    V <- V %*% sv$v
    d <- as.numeric(sv$d)

    if (!is.null(d_prev) && polish_tol > 0 && length(d_prev) == length(d)) {
      rel_change <- max(abs(d - d_prev) / pmax(abs(d_prev), 1e-12))
      if (is.finite(rel_change) && rel_change < polish_tol) {
        break
      }
    }
    d_prev <- d
  }

  list(u = U, v = V, d = d)
}

# Randomized low-pass GMD backend in pure R (fallback):
# 2 passes when n_power = 0, 2 + 2*n_power passes otherwise.
# `tol` is a relative singular-value cutoff (d_j <= tol * d_1 is dropped) and
# `jitter` is relative to the Gram scale in metric_orthonormalize().
gmd_randomized_r <- function(X, Q, R, k,
                             oversample = 20L,
                             n_power = 1L,
                             n_polish = 0L,
                             jitter = 1e-10,
                             tol = 1e-6,
                             polish_tol = 0,
                             seed = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  if (!is.numeric(k) || length(k) != 1 || k < 1 || k != floor(k)) {
    stop("k must be a positive integer.")
  }
  if (!is.numeric(oversample) || length(oversample) != 1 || oversample < 0 || oversample != floor(oversample)) {
    stop("oversample must be a non-negative integer.")
  }
  if (!is.numeric(n_power) || length(n_power) != 1 || n_power < 0 || n_power != floor(n_power)) {
    stop("n_power must be a non-negative integer.")
  }
  if (!is.numeric(n_polish) || length(n_polish) != 1 || n_polish < 0 || n_polish != floor(n_polish)) {
    stop("n_polish must be a non-negative integer.")
  }
  if (!is.numeric(polish_tol) || length(polish_tol) != 1 || polish_tol < 0) {
    stop("polish_tol must be a single non-negative number.")
  }

  k <- as.integer(min(k, n, p))
  ell <- as.integer(min(n, p, k + oversample))
  empty <- list(u = matrix(0.0, nrow = n, ncol = 0),
                v = matrix(0.0, nrow = p, ncol = 0),
                d = numeric(0), k = 0L)
  if (ell < 1) return(empty)

  if (!inherits(Q, "Matrix")) Q <- Matrix::Matrix(Q, sparse = FALSE)
  if (!inherits(R, "Matrix")) R <- Matrix::Matrix(R, sparse = FALSE)
  applyQ <- function(B) Q %*% B
  applyR <- function(B) R %*% B

  if (ell == n) {
    U0 <- metric_orthonormalize(diag(n), applyQ, jitter = jitter)
  } else {
    Omega <- if (ell == p) diag(p) else random_normal_matrix(p, ell, seed = seed)
    Y <- X %*% applyR(Omega)
    if (n_power > 0L) {
      for (it in seq_len(as.integer(n_power))) {
        Utmp <- metric_orthonormalize(Y, applyQ, jitter = jitter)
        Z <- Matrix::crossprod(X, applyQ(Utmp))
        Y <- X %*% applyR(Z)
      }
    }
    U0 <- metric_orthonormalize(Y, applyQ, jitter = jitter)
  }
  if (ncol(U0) == 0) return(empty)

  B <- Matrix::crossprod(X, applyQ(U0))
  RB <- applyR(B)
  G <- Matrix::crossprod(B, RB)
  G <- 0.5 * (G + t(G))

  eg <- eigen(as.matrix(G), symmetric = TRUE)
  k_use <- min(k, ncol(eg$vectors))
  if (k_use < 1) return(empty)

  S <- eg$vectors[, seq_len(k_use), drop = FALSE]
  d <- sqrt(pmax(eg$values[seq_len(k_use)], 0))

  U <- U0 %*% S
  V <- B %*% S
  nz <- .keep_components(d, tol)
  if (any(nz)) {
    V[, nz] <- sweep(V[, nz, drop = FALSE], 2, d[nz], "/")
  }
  if (any(!nz)) {
    V[, !nz] <- 0.0
  }

  if (n_polish > 0L) {
    pol <- gmd_randomized_polish(
      X = X,
      applyQ = applyQ,
      applyR = applyR,
      U = U,
      V = V,
      iters = as.integer(n_polish),
      jitter = jitter,
      polish_tol = polish_tol
    )
    U <- pol$u
    V <- pol$v
    d <- pol$d
  }

  keep <- .keep_components(d, tol)
  if (!any(keep)) return(empty)

  U <- U[, keep, drop = FALSE]
  V <- V[, keep, drop = FALSE]
  d <- as.numeric(d[keep])

  list(u = U, v = V, d = d, k = length(d))
}

gmd_randomized <- function(X, Q, R, k,
                           oversample = 20L,
                           n_power = 1L,
                           n_polish = 0L,
                           jitter = 1e-10,
                           tol = 1e-6,
                           polish_tol = 0,
                           seed = NULL,
                           use_cpp = TRUE) {
  if (isTRUE(use_cpp) && exists("gmd_randomized_cpp_dn", mode = "function")) {
    if (!inherits(Q, "Matrix")) Q <- Matrix::Matrix(Q, sparse = FALSE)
    if (!inherits(R, "Matrix")) R <- Matrix::Matrix(R, sparse = FALSE)

    seed_val <- if (is.null(seed)) 1234L else as.integer(seed)
    common <- list(X = X, k = as.integer(k), oversample = as.integer(oversample),
                   n_power = as.integer(n_power), n_polish = as.integer(n_polish),
                   jitter = jitter, tol = tol, polish_tol = polish_tol, seed = seed_val)
    cpp_res <- tryCatch({
      if (methods::is(Q, "sparseMatrix") && methods::is(R, "sparseMatrix")) {
        do.call(gmd_randomized_cpp_sp, c(list(Q = as_dgc(Q), R = as_dgc(R)), common))
      } else if (methods::is(Q, "sparseMatrix")) {
        do.call(gmd_randomized_cpp_qsp_rdn, c(list(Q = as_dgc(Q), R = as.matrix(R)), common))
      } else if (methods::is(R, "sparseMatrix")) {
        do.call(gmd_randomized_cpp_qdn_rsp, c(list(Q = as.matrix(Q), R = as_dgc(R)), common))
      } else {
        do.call(gmd_randomized_cpp_dn, c(list(Q = as.matrix(Q), R = as.matrix(R)), common))
      }
    }, error = function(e) {
      warning("C++ randomized backend failed, falling back to R implementation: ", e$message)
      NULL
    })

    if (!is.null(cpp_res)) {
      if (is.matrix(cpp_res$d)) cpp_res$d <- as.vector(cpp_res$d)
      cpp_res$k <- length(cpp_res$d)
      return(cpp_res)
    }
  }

  gmd_randomized_r(
    X = X,
    Q = Q,
    R = R,
    k = k,
    oversample = oversample,
    n_power = n_power,
    n_polish = n_polish,
    jitter = jitter,
    tol = tol,
    polish_tol = polish_tol,
    seed = seed
  )
}
