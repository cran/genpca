#' Internal metric and operator helpers
#'
#' These utilities are used across `genpls()` and `gplssvd_op()` to construct
#' metric-aware linear operators for GPLSSVD/PLS-SVD without duplicating logic.
#'
#' @details
#' - `.metric_operators(W, n_expected, remedy, name)` accepts a symmetric
#'   positive (semi)definite matrix or a diagonal specification and returns a
#'   list of closures that apply `W`, `W^{1/2}`, and `W^{-1/2}` to
#'   vectors/matrices. Identity and diagonal cases use fast paths; general PSD
#'   uses an eigen-based symmetric square root (pseudo-inverse on the null
#'   space). Input is validated; an indefinite weight is an error unless
#'   `remedy` asks for a repair, which is reported with a
#'   `genpca_metric_repaired` warning (see [repair_metric()]).
#'   The list also carries `kind` (one of `"identity"`, `"diag"`, `"general"`)
#'   and `mult_sqrt_right(x)` computing `x %*% W^{1/2}`, so callers can pick
#'   sparsity-preserving strategies.
#' - `.build_pls_operator(X, Y, MX, MY, WX, WY)` builds the implicit operator
#'   `v -> WX^{1/2} X^T MX^{1/2} MY^{1/2} Y WY^{1/2} v` and its adjoint. When the
#'   whitened blocks `Xe = MX^{1/2} X WX^{1/2}` and `Ye = MY^{1/2} Y WY^{1/2}`
#'   can be materialized without densifying sparse data, they are precomputed
#'   once so each matrix-vector product costs two multiplies instead of six.
#'
#' These helpers assume inputs are conformable and (when required) symmetric PSD;
#' callers are responsible for any validation or coercion.
#'
#' @keywords internal
#' @noRd

# Build metric operators for a symmetric (PS) matrix W:
# returns closures: mult (W %*% x), mult_sqrt (W^{1/2} x), mult_invsqrt (W^{-1/2} x),
# mult_sqrt_right (x %*% W^{1/2}), plus a `kind` tag.
.metric_operators <- function(W, n_expected = NULL,
                              remedy = c("error", "ridge", "clip", "identity"),
                              name = "weight matrix") {
  remedy <- match.arg(remedy)
  rtol <- .metric_rtol_default()
  repair <- function(M) {
    if (remedy == "error") {
      stop(name, " must be symmetric positive semi-definite", call. = FALSE)
    }
    B <- repair_metric(M, method = remedy, rtol = rtol, name = name)
    if (isTRUE(attr(B, "repair_report")$changed)) .warn_metric_repaired(attr(B, "repair_report"))
    attr(B, "repair_report") <- NULL
    B
  }
  # Identity
  if (is.null(W)) {
    id <- function(x) x
    return(list(
      kind = "identity",
      mult = id,
      mult_sqrt = id,
      mult_invsqrt = id,
      mult_sqrt_right = id
    ))
  }

  # Diagonal supplied as a numeric vector
  if (is.numeric(W) && is.null(dim(W))) {
    if (!is.null(n_expected) && length(W) != n_expected) {
      stop("numeric weight vector has length ", length(W),
           " but length ", n_expected, " is required")
    }
    d <- as.numeric(W)
  } else if ((is.matrix(W) || inherits(W, "Matrix")) && Matrix::isDiagonal(W)) {
    if (!is.null(n_expected) && nrow(W) != n_expected) {
      stop("weight matrix is ", nrow(W), " x ", ncol(W),
           " but dimension ", n_expected, " is required")
    }
    d <- as.numeric(Matrix::diag(W))
  } else {
    d <- NULL
  }

  # Diagonal fast path: `v * x` recycles column-wise, i.e. row scaling for
  # matrices (dense or sparse) and plain elementwise for vectors.
  if (!is.null(d)) {
    if (any(!is.finite(d))) stop("Diagonal elements of ", name, " must be finite", call. = FALSE)
    if ((remedy == "clip" && any(d < 0)) ||
        (length(d) && max(abs(d)) > 0 && any(d < -rtol * max(abs(d))))) {
      d <- as.numeric(Matrix::diag(repair(Matrix::Diagonal(length(d), x = d))))
    }
    d <- .clamp_weights(d, name = name)
    ds <- sqrt(d)
    invds <- ifelse(d > 0, 1 / ds, 0)
    return(list(
      kind = "diag",
      mult = function(x) d * x,
      mult_sqrt = function(x) ds * x,
      mult_invsqrt = function(x) invds * x,
      mult_sqrt_right = function(x) {
        if (inherits(x, "Matrix")) x %*% Matrix::Diagonal(x = ds)
        else x * rep(ds, each = nrow(x))
      }
    ))
  }

  # General symmetric PSD -> symmetric sqrt via eigen
  if (!is.null(n_expected) && (nrow(W) != n_expected || ncol(W) != n_expected)) {
    stop("weight matrix is ", nrow(W), " x ", ncol(W),
         " but dimension ", n_expected, " is required")
  }
  # Validate, never repair: asymmetry beyond roundoff and indefiniteness
  # are input errors here (callers that offer a remedy apply it first).
  W <- symmetrize_or_stop(W, name = name)
  if (remedy == "clip" || !is_psd(W, rtol = rtol)) W <- repair(W)
  Wd <- as.matrix(W)
  es <- eigen(Wd, symmetric = TRUE)
  lam <- pmax(es$values, 0)
  # Eigenvalues within metric_rtol of the largest are the null space: W^{1/2}
  # and the pseudo-inverse W^{-1/2} both vanish there.
  cut <- .metric_rtol_default() * max(lam, 0)
  # dgeMatrix products call BLAS directly (no base-R NaN scan per product)
  Q <- Matrix::Matrix(es$vectors, sparse = FALSE)
  sqrt_lam <- ifelse(lam > cut, sqrt(lam), 0)
  invsqrt_lam <- ifelse(lam > cut, 1 / sqrt(lam), 0)
  apply_factored <- function(x, scal) {
    alpha <- scal * Matrix::crossprod(Q, x) # row scaling, no k x k diag matmul
    Q %*% alpha
  }
  list(
    kind = "general",
    mult = function(x) W %*% x,
    mult_sqrt = function(x) apply_factored(x, sqrt_lam),
    mult_invsqrt = function(x) apply_factored(x, invsqrt_lam),
    # W^{1/2} is symmetric, so x %*% W^{1/2} = t(W^{1/2} t(x))
    mult_sqrt_right = function(x) t(apply_factored(t(x), sqrt_lam))
  )
}

#' Build PLS operator closures (internal)
#'
#' @param X Matrix `n x p`
#' @param Y Matrix `n x q`
#' @param MX,MY,WX,WY Metric operator lists from `.metric_operators()`
#' @return A list with `S_mv`, `ST_mv`, `dims = c(p, q)`, `materialized`, and
#'   (when materialized) the whitened blocks `Xe`, `Ye`.
#' @keywords internal
#' @noRd
#'
# Build PLS operator closures S_mv and ST_mv for given data and metric operators
.build_pls_operator <- function(X, Y, MX, MY, WX, WY) {
  # X: n x p, Y: n x q; MX,MY,WX,WY are metric operator lists (mult/mult_sqrt/...)
  p <- ncol(X)
  q <- ncol(Y)

  # Materializing Xe = MX^{1/2} X WX^{1/2} is safe unless X is sparse and a
  # metric is dense-general (which would densify X). Scaling metrics preserve
  # sparsity, and for dense X the whitened copy has the same footprint as X.
  can_mat <- function(A, Mrow, Wcol) {
    !inherits(A, "sparseMatrix") ||
      (Mrow$kind != "general" && Wcol$kind != "general")
  }

  if (can_mat(X, MX, WX) && can_mat(Y, MY, WY)) {
    # Store dense blocks as dgeMatrix: Matrix products call BLAS directly,
    # skipping base R's per-product NaN scan of the full matrix.
    as_dge <- function(A) if (inherits(A, "Matrix")) A else Matrix::Matrix(A, sparse = FALSE)
    Xe <- as_dge(WX$mult_sqrt_right(MX$mult_sqrt(X)))
    Ye <- as_dge(WY$mult_sqrt_right(MY$mult_sqrt(Y)))
    S_mv  <- function(v, args = NULL) as.matrix(Matrix::crossprod(Xe, Ye %*% v))
    ST_mv <- function(u, args = NULL) as.matrix(Matrix::crossprod(Ye, Xe %*% u))
    return(list(S_mv = S_mv, ST_mv = ST_mv, dims = c(p, q),
                materialized = TRUE, Xe = Xe, Ye = Ye))
  }

  # Lazy path: apply the whitening chain per product. When both blocks share
  # the same row metric, MX^{1/2} MY^{1/2} = MX, so the middle collapses to a
  # single multiply (and avoids compounding square-root error).
  mid <- if (identical(MX, MY)) MX$mult else function(x) MX$mult_sqrt(MY$mult_sqrt(x))
  mid_t <- if (identical(MX, MY)) MX$mult else function(x) MY$mult_sqrt(MX$mult_sqrt(x))
  S_mv <- function(v, args = NULL) {
    t1 <- Y %*% WY$mult_sqrt(v)
    t2 <- Matrix::crossprod(X, mid(t1))
    as.matrix(WX$mult_sqrt(t2))
  }
  ST_mv <- function(u, args = NULL) {
    s1 <- X %*% WX$mult_sqrt(u)
    s2 <- Matrix::crossprod(Y, mid_t(s1))
    as.matrix(WY$mult_sqrt(s2))
  }
  list(S_mv = S_mv, ST_mv = ST_mv, dims = c(p, q), materialized = FALSE)
}
