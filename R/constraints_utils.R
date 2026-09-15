#' @title Utilities for constraints
#' @description Helpers to validate, symmetrize and (only when asked) repair
#' constraint matrices. Two relative tolerances are used throughout the
#' package: `metric_rtol` (default `sqrt(.Machine$double.eps)`) decides
#' whether a metric is positive (semi)definite and which of its eigenvalues
#' count as zero, and `rank_rtol` (see [genpca()]) decides which components
#' are kept. Both are relative to the scale of the matrix, so every decision
#' is invariant to rescaling the input.
#' @name constraints_utils
#' @keywords internal
NULL

#' Default relative tolerance for metric validation
#' @keywords internal
.metric_rtol_default <- function() sqrt(.Machine$double.eps)

# Scale used for relative tolerances on a symmetric matrix: max |diag|.
# Returns 1 for the zero matrix and NA when the diagonal is not finite.
.metric_scale <- function(A) {
  d <- abs(as.numeric(Matrix::diag(A)))
  s <- if (length(d)) max(d) else 0
  if (!is.finite(s)) return(NA_real_)
  if (s == 0) {
    # A PSD matrix has its largest entry on the diagonal, so a zero diagonal
    # with non-zero off-diagonals is indefinite at the off-diagonal scale.
    off <- if (methods::is(A, "sparseMatrix")) {
      if (length(A@x)) max(abs(A@x)) else 0
    } else {
      max(abs(as.matrix(A)), 0)
    }
    if (!is.finite(off)) return(NA_real_)
    s <- off
  }
  if (s == 0) 1 else s
}

.as_metric_matrix <- function(A) {
  if (inherits(A, "Matrix")) A else Matrix::Matrix(A, sparse = FALSE)
}

# TRUE iff A + shift * I admits a Cholesky factorization. The factorization
# must FAIL on indefinite input: dense matrices go through base::chol
# (dpotrf) because Matrix::Cholesky() on dense input uses pivoted dpstrf,
# which succeeds with a warning on indefinite matrices.
.chol_probe <- function(A, shift = 0) {
  probe <- if (shift != 0) A + Matrix::Diagonal(nrow(A), x = shift) else A
  if (methods::is(probe, "sparseMatrix")) {
    suppressWarnings(tryCatch({
      Matrix::Cholesky(probe, LDL = FALSE, Imult = 0, super = TRUE)
      TRUE
    }, error = function(e) FALSE))
  } else {
    tryCatch({
      chol(as.matrix(probe))
      TRUE
    }, error = function(e) FALSE)
  }
}

#' @title Test positive semi-definiteness (relative tolerance)
#' @description `is_psd()` is TRUE when `A` is symmetric and every
#' eigenvalue exceeds `-rtol * max(abs(diag(A)))`; `is_pd()` is TRUE when
#' every eigenvalue exceeds `+rtol * max(abs(diag(A)))`. Both are shifted
#' Cholesky probes, so large sparse matrices never need an
#' eigendecomposition. `is_spd()` is a deprecated alias of `is_psd()` kept
#' for internal callers (its `tol` is the relative tolerance).
#' @param A numeric matrix or Matrix::Matrix
#' @param rtol relative tolerance (default `sqrt(.Machine$double.eps)`)
#' @return logical
#' @keywords internal
is_psd <- function(A, rtol = .metric_rtol_default()) {
  A <- .as_metric_matrix(A)
  if (nrow(A) != ncol(A) || !Matrix::isSymmetric(A)) return(FALSE)
  s <- .metric_scale(A)
  if (is.na(s)) return(FALSE)
  .chol_probe(A, rtol * s)
}

#' @rdname is_psd
#' @keywords internal
is_pd <- function(A, rtol = .metric_rtol_default()) {
  A <- .as_metric_matrix(A)
  if (nrow(A) != ncol(A) || !Matrix::isSymmetric(A)) return(FALSE)
  s <- .metric_scale(A)
  if (is.na(s)) return(FALSE)
  .chol_probe(A, -rtol * s)
}

#' @rdname is_psd
#' @param tol relative tolerance (deprecated name; same as `rtol`)
#' @keywords internal
is_spd <- function(A, tol = .metric_rtol_default()) {
  is_psd(A, rtol = tol)
}

#' @title Symmetrize a nearly symmetric matrix or stop
#' @description Measures `||A - A'||_F / ||A||_F`. Below `rtol` the two
#' triangles are averaged and the result is marked symmetric; above it the
#' function stops. Asymmetry is an input error, never something a PSD remedy
#' repairs.
#' @param A square numeric matrix or Matrix::Matrix
#' @param rtol relative asymmetry allowed (default 1e-10)
#' @param name label used in the error message
#' @return a symmetric Matrix (dense or sparse as supplied)
#' @keywords internal
symmetrize_or_stop <- function(A, rtol = 1e-10, name = "A") {
  A <- .as_metric_matrix(A)
  if (nrow(A) != ncol(A)) {
    stop("Matrix ", name, " must be square", call. = FALSE)
  }
  if (methods::is(A, "symmetricMatrix") || methods::is(A, "diagonalMatrix")) {
    return(A)
  }
  D <- A - Matrix::t(A)
  nD <- Matrix::norm(D, "F")
  if (nD == 0) return(Matrix::forceSymmetric(A))
  nA <- Matrix::norm(A, "F")
  rel <- nD / max(nA, .Machine$double.xmin)
  if (is.finite(rel) && rel <= rtol) {
    return(Matrix::forceSymmetric((A + Matrix::t(A)) / 2))
  }
  stop("Matrix ", name, " must be symmetric (relative asymmetry ",
       signif(rel, 3), " exceeds ", rtol, ")", call. = FALSE)
}

#' @title Validate diagonal weights
#' @description Weights must be finite and non-negative. Entries in
#' `[-rtol * max(abs(w)), 0)` are set to exactly zero with a message;
#' anything below is an error. Every backend sees the same cleaned vector.
#' @keywords internal
.clamp_weights <- function(w, rtol = .metric_rtol_default(), name = "A") {
  w <- as.numeric(w)
  if (any(!is.finite(w))) {
    stop("Diagonal elements of ", name, " must be finite", call. = FALSE)
  }
  if (!length(w)) return(w)
  s <- max(abs(w))
  if (s == 0) return(w)
  if (any(w < -rtol * s)) {
    stop("Diagonal elements of ", name, " must be non-negative (smallest = ",
         signif(min(w), 3), ")", call. = FALSE)
  }
  if (any(w < 0)) {
    message("Setting ", sum(w < 0), " tiny negative weight(s) in ", name,
            " to zero (all within ", signif(rtol, 3), " * max weight).")
    w[w < 0] <- 0
  }
  w
}

#' @title Coerce to general CSC sparse matrix (dgCMatrix)
#' @description Replacement for the deprecated direct `as(., "dgCMatrix")`
#' coercion from symmetric/triangular/diagonal Matrix classes.
#' @param A a matrix or Matrix
#' @return a dgCMatrix
#' @keywords internal
as_dgc <- function(A) {
  if (inherits(A, "dgCMatrix")) return(A)
  if (!inherits(A, "Matrix")) A <- Matrix::Matrix(A, sparse = TRUE)
  methods::as(methods::as(methods::as(A, "dMatrix"), "generalMatrix"), "CsparseMatrix")
}

#' @title Coerce dense symmetric Matrix classes to general dense (dgeMatrix)
#' @description Replacement for the deprecated direct `as(., "dgeMatrix")`
#' coercion from dsyMatrix/dpoMatrix.
#' @param A a dense Matrix
#' @return a dgeMatrix
#' @keywords internal
as_dge <- function(A) {
  if (inherits(A, "dgeMatrix")) return(A)
  methods::as(methods::as(A, "generalMatrix"), "unpackedMatrix")
}

#' @title Clip a symmetric matrix to the PSD cone
#' @description Spectral clip: eigen-decompose and set negative eigenvalues to
#' zero. Unlike [ensure_spd()] (a diagonal ridge shift), this preserves the
#' non-negative part of the spectrum exactly. The output has no negative
#' eigenvalue beyond reconstruction roundoff: the only fast path is an exact
#' (unshifted) Cholesky success, which proves positive definiteness. Requires a dense
#' eigendecomposition, so large sparse matrices are refused.
#' @param M numeric matrix or Matrix::Matrix
#' @param tol unused; kept for call compatibility
#' @param dense_maxn refuse sparse input larger than this (clip densifies)
#' @param name label used in error messages
#' @return a symmetric Matrix, PSD
#' @keywords internal
clip_psd <- function(M, tol = NULL, dense_maxn = 2000L, name = "M") {
  M <- symmetrize_or_stop(M, name = name)
  n <- nrow(M)
  if (Matrix::isDiagonal(M)) {
    return(Matrix::Diagonal(n, x = pmax(as.numeric(Matrix::diag(M)), 0)))
  }
  if (is_pd(M, rtol = 0)) return(M)
  if (methods::is(M, "sparseMatrix") && n > dense_maxn) {
    stop("constraints_remedy = 'clip' needs a dense eigendecomposition, but ",
         "the matrix is sparse with ", n, " rows (> ", dense_maxn,
         "). Use constraints_remedy = 'ridge' instead.")
  }
  ee <- eigen(as.matrix(M), symmetric = TRUE)
  vals <- pmax(ee$values, 0)
  Mc <- ee$vectors %*% (vals * t(ee$vectors))
  Matrix::forceSymmetric(Matrix::Matrix((Mc + t(Mc)) / 2, sparse = FALSE))
}

#' @title Ensure SPD (sparse-friendly)
#' @description Force a symmetric matrix to be symmetric positive definite:
#' the result satisfies `is_pd(., rtol = tol)`. Already-PD input is returned
#' unchanged; otherwise a Gershgorin-based diagonal shift is applied, with a
#' `Matrix::nearPD()` fallback for small dense matrices and an escalating
#' jitter as a last resort.
#' @param M numeric matrix or Matrix::Matrix
#' @param tol relative positive-definiteness margin (default 1e-6)
#' @param nearpd_maxn only use nearPD when n <= nearpd_maxn and matrix is dense
#' @param name label used in error messages
#' @return a Matrix object (sparse stays sparse when possible)
#' @keywords internal
#' @importFrom Matrix forceSymmetric Diagonal rowSums Cholesky nearPD
#' @importFrom methods as is
ensure_spd <- function(M, tol = 1e-6, nearpd_maxn = 2000L, name = "M") {
  M <- .as_metric_matrix(M)
  xvals <- if (methods::is(M, "sparseMatrix")) M@x else as.numeric(M)
  if (length(xvals) && !all(is.finite(xvals))) {
    stop("ensure_spd: matrix contains non-finite values.")
  }
  M <- symmetrize_or_stop(M, name = name)
  n <- nrow(M)
  scale <- .metric_scale(M)

  if (is_pd(M, rtol = tol)) return(M)

  # Gershgorin circle theorem: eigenvalues lie in circles centred at a_ii
  # with radius sum_{j != i} |a_ij|. A diagonal shift that makes every
  # a_ii - r_i positive guarantees positive definiteness.
  d <- Matrix::diag(M)
  rs <- Matrix::rowSums(abs(M)) - abs(d)
  min_margin <- suppressWarnings(min(d - rs))
  if (!is.finite(min_margin)) min_margin <- -1
  if (min_margin <= tol * scale) {
    shift <- (-min_margin) + max(tol * scale, abs(min_margin) * 0.1)
    M <- M + Matrix::Diagonal(n, x = shift)
  }
  if (is_pd(M, rtol = tol)) return(M)

  if (!methods::is(M, "sparseMatrix") && n <= nearpd_maxn) {
    Mnp <- Matrix::nearPD(as.matrix(M), corr = FALSE, keepDiag = TRUE)$mat
    Mnp <- Matrix::forceSymmetric(Matrix::Matrix(as.matrix(Mnp), sparse = FALSE))
    if (is_pd(Mnp, rtol = tol)) return(Mnp)
    M <- Mnp
  }

  jitter <- max(tol * scale, 1e-10 * scale)
  for (i in seq_len(10)) {
    M2 <- M + Matrix::Diagonal(n, x = jitter)
    if (is_pd(M2, rtol = tol)) return(M2)
    jitter <- jitter * 100
  }

  stop("ensure_spd: Unable to make matrix SPD; consider revising constraint.")
}
