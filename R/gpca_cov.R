#' Generalized PCA on a covariance matrix (GMD form)
#'
#' Performs Generalized PCA directly on a pre-computed covariance matrix
#' `C = X'MX` with a single variable-side metric `R`, following Allen et
#' al.'s GMD: the eigendecomposition of \eqn{R^{1/2} C R^{1/2}} mapped back
#' with \eqn{V = R^{-1/2} Z}, so that \eqn{V'RV = I}. With `C = X'MX` and
#' `R = A` this matches \code{\link{genpca}(X, M = M, A = A)} exactly. This is
#' useful when you already have `C` or when `X` is too large to store but `C`
#' is manageable.
#'
#' The generalized eigenproblem \eqn{C v = \lambda R v} is a different
#' estimator (it maximises \eqn{v'Cv} subject to \eqn{v'Rv = 1}, which is
#' generally gives different components from the GMD) and lives in its own
#' function, \code{\link{geigen_cov}}. `method = "geigen"` is accepted here
#' for one release and forwards to it with a deprecation warning.
#'
#' @param C A p x p symmetric positive semi-definite covariance matrix,
#'   typically `C = X'MX`. Asymmetry beyond roundoff and indefiniteness
#'   beyond `metric_rtol` are errors.
#' @param R Variable-side constraint/metric. Can be:
#'   \itemize{
#'     \item NULL: identity matrix (standard PCA on C)
#'     \item a numeric vector of length p: diagonal weights (must be non-negative)
#'     \item a p x p symmetric PSD matrix: general metric/smoothing/structure penalties
#'   }
#' @param ncomp Number of components to return. Default is all positive eigenvalues.
#' @param method Deprecated. `"gmd"` (default) is this function; `"geigen"`
#'   forwards to \code{\link{geigen_cov}} with a warning.
#' @param constraints_remedy Deprecated here (GMD requires PSD input and stops
#'   otherwise); forwarded to \code{\link{geigen_cov}} when
#'   `method = "geigen"`.
#' @param rank_rtol Relative cutoff for component acceptance on the
#'   singular-value scale (components with `d_j <= rank_rtol * d_1`
#'   are dropped). Default 1e-6.
#' @param metric_rtol Relative tolerance for validating \code{C} and
#'   \code{R} and for detecting the numerical null space in an
#'   eigendecomposition of a general \code{R}. Every strictly positive
#'   diagonal weight is retained without a rank approximation. Default
#'   \code{sqrt(.Machine$double.eps)}.
#' @param tol Deprecated; use \code{rank_rtol} and \code{metric_rtol}.
#' @param verbose Logical. If TRUE, print progress messages. Default FALSE.
#'
#' @return A plain list (\strong{not} a \pkg{multivarious}
#'   \code{bi_projector}) with components:
#'   \describe{
#'     \item{v}{p x k matrix of loadings (R-orthonormal eigenvectors)}
#'     \item{d}{Singular values (square root of eigenvalues lambda)}
#'     \item{lambda}{Eigenvalues (variances under the R-metric)}
#'     \item{k}{Number of components returned}
#'     \item{propv}{Proportion of variance explained by each component
#'       (total variance is \eqn{\mathrm{tr}(CR)}, Allen et al. Corollary 5)}
#'     \item{cumv}{Cumulative proportion of variance explained}
#'     \item{R_rank}{Rank of the constraint matrix R}
#'     \item{method}{`"gmd"`}
#'   }
#'   Because this is a plain list rather than a \code{bi_projector}, the
#'   \code{multivarious} generics \code{scores()}, \code{components()}, and
#'   \code{reconstruct()} do not apply to it; index \code{$v}/\code{$d}
#'   directly, or use \code{\link{genpca}} when you need the full projector
#'   interface on a data matrix rather than a pre-computed covariance matrix.
#'
#' @examples
#' # Standard PCA on a covariance (no constraint)
#' C <- cov(scale(iris[,1:4], center=TRUE, scale=FALSE))
#' fit0 <- genpca_cov(C, R=NULL, ncomp=3)
#' print(fit0$d[1:3])       # first 3 singular values
#' print(fit0$propv[1:3])   # variance explained by first 3 components
#'
#' # Equivalence with genpca()
#' set.seed(123)
#' X <- matrix(rnorm(50 * 10), 50, 10)
#' M_diag <- runif(50, 0.5, 1.5)  # row weights
#' A_diag <- runif(10, 0.5, 2)    # column weights
#' fit_gpca <- genpca(X, M = M_diag, A = A_diag, ncomp = 5,
#'                    preproc = multivarious::pass())
#' C <- crossprod(X, diag(M_diag) %*% X)  # C = X'MX
#' fit_cov <- genpca_cov(C, R = A_diag, ncomp = 5)
#' all.equal(fit_gpca$sdev, fit_cov$d, tolerance = 1e-10)
#'
#' # Variable weights via a diagonal metric (iris covariance, 4 variables)
#' C_iris <- cov(scale(iris[,1:4], center=TRUE, scale=FALSE))
#' w <- c(1, 1, 0.5, 2)
#' fitW <- genpca_cov(C_iris, R = w, ncomp=3)
#' print(fitW$d[1:3])
#'
#' @seealso \code{\link{geigen_cov}} for the generalized eigenproblem
#'   \eqn{C v = \lambda R v}, \code{\link{genpca}} for the two-sided GPCA on
#'   data matrices, \code{\link{genpls}} for generalized partial least squares
#'
#' @references
#' Allen, G. I., Grosenick, L., & Taylor, J. (2014).
#' A Generalized Least-Squares Matrix Decomposition.
#' Journal of the American Statistical Association, 109(505), 145-159.
#'
#' @export
#' @importFrom Matrix Matrix isSymmetric forceSymmetric Diagonal t diag crossprod
genpca_cov <- function(C, R = NULL, ncomp = NULL,
                       method = c("gmd", "geigen"),
                       constraints_remedy = c("error", "ridge", "clip", "identity"),
                       rank_rtol = 1e-6, metric_rtol = .metric_rtol_default(),
                       tol = NULL, verbose = FALSE) {

  method <- match.arg(method)
  if (!is.null(tol)) {
    warning("`tol` is deprecated in genpca_cov(); use `rank_rtol` and `metric_rtol`.", call. = FALSE)
  }
  remedy_supplied <- !missing(constraints_remedy)
  constraints_remedy <- match.arg(constraints_remedy)

  if (method == "geigen") {
    warning("genpca_cov(method = \"geigen\") is deprecated; call geigen_cov() directly.", call. = FALSE)
    return(geigen_cov(C, R, ncomp, constraints_remedy = constraints_remedy,
                      rank_rtol = rank_rtol, metric_rtol = metric_rtol, verbose = verbose))
  }
  if (remedy_supplied && constraints_remedy != "error") {
    warning("`constraints_remedy` is ignored by genpca_cov(): the GMD form requires PSD input. ",
            "Repair the metric explicitly with repair_metric(), or use geigen_cov().", call. = FALSE)
  }
  genpca_cov_gmd(C, R, ncomp, rank_rtol, metric_rtol, verbose)
}

#' Generalized eigenproblem on a covariance matrix
#'
#' Maximises \eqn{v'Cv} subject to \eqn{v'Rv = 1} and the additional
#' constraint that \eqn{v} lies in the retained range of `R`. Successive
#' components are \eqn{R}-orthogonal. If \eqn{P} is the orthogonal projector
#' onto that range, the returned vectors satisfy
#' \eqn{P C v = \lambda R v} and \eqn{V'RV = I}. For a full-rank `R` this is
#' the usual equation \eqn{C v = \lambda R v}. It also holds for a singular
#' `R` when `C` maps its retained range into itself. Otherwise the component
#' of \eqn{C v} outside the retained range need not vanish.
#'
#' This is a different estimator from the GMD of \code{\link{genpca_cov}},
#' which uses \eqn{R^{1/2} C R^{1/2}}. If `C` and `R` commute, their common
#' eigenvectors can be ordered differently: the GMD weights variances by
#' metric eigenvalues, whereas this estimator divides by them. With
#' `R = c * I`, the directions and their ordering agree, but the eigenvalue
#' scales differ unless `c = 1`.
#'
#' `C` is validated for symmetry but may be indefinite (the generalized
#' eigenproblem is still defined); a warning is issued when its minimum
#' eigenvalue is below `-metric_rtol * scale`. `R` must be positive
#' semi-definite; an indefinite `R` is subject to `constraints_remedy`.
#'
#' @inheritParams genpca_cov
#' @param C A p x p symmetric matrix. Asymmetry beyond roundoff is an error;
#'   an indefinite `C` is allowed (the problem is still defined) and only
#'   produces a warning.
#' @param constraints_remedy What to do with an indefinite `R`: `"error"`
#'   (default), `"ridge"`, `"clip"` or `"identity"`; a repair emits a
#'   `genpca_metric_repaired` warning. See \code{\link{genpca}}.
#' @return A plain list with the same components as \code{\link{genpca_cov}}
#'   (`v`, `d`, `lambda`, `k`, `propv`, `cumv`, `R_rank`) and
#'   `method = "geigen"`. `propv` is relative to
#'   \eqn{\mathrm{tr}(R^{-1/2} C R^{-1/2})} on the range of `R`.
#' @examples
#' C <- cov(scale(iris[,1:4], center=TRUE, scale=FALSE))
#' w <- c(1, 1, 0.5, 2)
#' fit_gmd <- genpca_cov(C, R = w, ncomp = 2)
#' fit_geigen <- geigen_cov(C, R = w, ncomp = 2)
#' # different estimators: the singular values generally differ
#' rbind(gmd = fit_gmd$d, geigen = fit_geigen$d)
#'
#' # With singular R, the equation is projected onto its retained range
#' C <- matrix(c(2, 1, 1, 2), 2)
#' R <- diag(c(1, 0))
#' fit <- geigen_cov(C, R, ncomp = 1)
#' P <- diag(c(1, 0))
#' P %*% C %*% fit$v - (R %*% fit$v) * fit$lambda
#' @seealso \code{\link{genpca_cov}}
#' @export
geigen_cov <- function(C, R = NULL, ncomp = NULL,
                       constraints_remedy = c("error", "ridge", "clip", "identity"),
                       rank_rtol = 1e-6, metric_rtol = .metric_rtol_default(),
                       verbose = FALSE) {
  constraints_remedy <- match.arg(constraints_remedy)
  genpca_cov_geigen(C, R, ncomp, constraints_remedy, rank_rtol, metric_rtol, verbose)
}

#' GMD-based covariance GPCA (internal)
#'
#' Implements Allen et al.'s GMD approach for covariance matrices.
#' Computes eigendecomposition of \eqn{R^{1/2} C R^{1/2}} and maps back.
#'
#' @keywords internal
#' @importFrom Matrix Matrix isSymmetric forceSymmetric Diagonal t diag
genpca_cov_gmd <- function(C, R = NULL, ncomp = NULL, rank_rtol = 1e-6,
                           metric_rtol = .metric_rtol_default(), verbose = FALSE) {

  # Basic checks & normalization. GMD models C = X' M X, so C must be PSD.
  stopifnot(is.matrix(C) || inherits(C, "Matrix"))
  p <- nrow(C)
  stopifnot(p == ncol(C))
  C <- symmetrize_or_stop(C, name = "C")
  if (!is_psd(C, rtol = metric_rtol)) {
    stop("C must be symmetric positive semi-definite", call. = FALSE)
  }

  # Column operator R (vector of weights, NULL=I, or PSD matrix)
  if (is.null(R)) {
    R <- Matrix::Diagonal(p)
  } else if (is.vector(R)) {
    stopifnot(length(R) == p)
    R <- Matrix::Diagonal(p, x = .clamp_weights(R, metric_rtol, "R"))
  } else {
    R <- symmetrize_or_stop(R, name = "R")
    if (!Matrix::isDiagonal(R) && !is_psd(R, rtol = metric_rtol)) {
      stop("R must be symmetric positive semi-definite", call. = FALSE)
    }
  }

  if (verbose) message("Computing eigen factorization of R...")

  # Factorization of R to build R^{1/2} and R^{-1/2} on range(R).
  # Diagonal R (including the weight-vector case) needs no eigendecomposition.
  if (Matrix::isDiagonal(R)) {
    r_diag <- .clamp_weights(as.numeric(Matrix::diag(R)), metric_rtol, "R")
    keep <- which(r_diag > 0)
    if (length(keep) == 0L) stop("R is (numerically) zero.")
    s_all <- sqrt(r_diag)
    Rsqrt <- Matrix::Diagonal(p, x = s_all)                                  # R^{1/2}
    Rmhalf <- Matrix::Diagonal(p, x = ifelse(r_diag > 0, 1 / s_all, 0))  # R^{-1/2} on range(R)
  } else {
    Re <- eigen(as.matrix(R), symmetric = TRUE)
    vals <- pmax(Re$values, 0)
    keep <- which(vals > metric_rtol * max(vals, 0))
    if (length(keep) == 0L) stop("R is (numerically) zero.")

    U <- Re$vectors[, keep, drop = FALSE]
    s <- sqrt(vals[keep])
    Rsqrt <- U %*% (s * t(U))                 # R^{1/2}
    Rmhalf <- U %*% ((1 / s) * t(U))           # R^{-1/2} on range(R)
  }

  if (verbose) message("Computing R^{1/2} C R^{1/2}...")

  # Core step per Allen: eigen of B = R^{1/2} C R^{1/2}
  B <- as.matrix(Rsqrt %*% C %*% Rsqrt)
  B <- 0.5 * (B + t(B))  # Ensure symmetry

  if (verbose) message("Computing eigendecomposition...")

  # Iterative top-k solver when few components are requested from a large B;
  # "LA" (largest algebraic) so tiny negative eigenvalues cannot be selected.
  if (!is.null(ncomp) && ncomp >= 1L && p > 100L && ncomp < (p - 1L)) {
    Ee <- .top_eigs_sym(B, ncomp, "LA", tol = 1e-10)
  } else {
    Ee <- eigen(B, symmetric = TRUE)
  }
  lam_all <- pmax(Ee$values, 0)

  # Keep components above the relative rank cutoff (on d^2 = lambda)
  pos <- which(lam_all > 0 & lam_all > rank_rtol^2 * max(lam_all, 0))
  if (length(pos) == 0L) stop("No positive eigenvalues in R^{1/2} C R^{1/2}.")
  if (is.null(ncomp)) ncomp <- length(pos)
  ncomp <- min(ncomp, length(pos))

  lam <- lam_all[pos][1:ncomp]
  Z <- Ee$vectors[, pos, drop = FALSE][, 1:ncomp, drop = FALSE]

  if (verbose) message("Mapping back: V = R^{-1/2} Z...")

  # Map back: V = R^{-1/2} Z  (so that V' R V = I)
  V <- Rmhalf %*% Z

  # Variance accounting in Allen's GPCA: total = tr(C R) = sum(C * R) (both symmetric)
  total <- sum(C * R)
  propv <- as.numeric(lam / ifelse(total > 0, total, 1))
  cumv <- cumsum(propv)

  list(
    v      = Matrix::Matrix(V, sparse = FALSE),  # p x k, R-orthonormal
    d      = sqrt(lam),                          # GMD values
    lambda = lam,                                # = d^2
    k      = ncomp,
    propv  = propv,
    cumv   = cumv,
    R_rank = length(keep),
    method = "gmd"
  )
}

#' Generalized eigenvalue-based covariance GPCA (internal)
#'
#' Solves the generalized eigenproblem projected onto the retained range of R.
#' This is the original implementation that was in gpca.R.
#'
#' @keywords internal
#' @importFrom Matrix Matrix isSymmetric forceSymmetric Diagonal t diag
genpca_cov_geigen <- function(C, R = NULL, ncomp = NULL,
                              constraints_remedy = c("error", "ridge", "clip", "identity"),
                              rank_rtol = 1e-6, metric_rtol = .metric_rtol_default(),
                              verbose = FALSE) {

  constraints_remedy <- match.arg(constraints_remedy)

  # --- Basic checks & normalization
  stopifnot(is.matrix(C) || inherits(C, "Matrix"))
  p <- nrow(C)
  stopifnot(p == ncol(C))
  C <- symmetrize_or_stop(C, name = "C")
  c_scale <- .metric_scale(C)

  # The generalized eigenproblem is defined for indefinite C; warn, do not stop.
  eigC_min <- tryCatch(
    if (p <= 800) min(eigen(as.matrix(C), symmetric = TRUE, only.values = TRUE)$values)
    else .top_eigs_sym(C, 1, "SA")$values,
    error = function(e) NA_real_
  )
  if (is.na(eigC_min)) {
    if (verbose) warning("Could not verify PSD of C; proceeding.")
  } else if (eigC_min < -metric_rtol * c_scale) {
    warning("C appears non-PSD (min eig: ", signif(eigC_min, 4), "). Proceeding but results may be unstable.")
  }

  # --- Prepare metric R: same validation and (reported) repair as genpca()
  R <- .prep_one_metric(R, p, "R", metric_rtol, constraints_remedy, verbose)

  if (verbose) message("Solving generalized eigenproblem C v = lambda R v...")

  # --- Solve C v = lambda R v in the range of R (handles semidefinite R)
  # Eigen-decompose R = U diag(gamma) U^T, keep gamma > tol.
  # Diagonal R needs no eigendecomposition: U is a sparse column selection.
  if (Matrix::isDiagonal(R)) {
    gamma <- pmax(as.numeric(Matrix::diag(R)), 0)
    keep <- which(gamma > 0)
    if (length(keep) == 0L) stop("R is numerically zero; no components can be extracted.")
    U <- Matrix::sparseMatrix(i = keep, j = seq_along(keep), x = 1,
                              dims = c(p, length(keep)))
    gam_sqrt_inv <- 1 / sqrt(gamma[keep])
  } else {
    ER <- eigen(as.matrix(R), symmetric = TRUE)
    gamma <- pmax(ER$values, 0)
    keep <- which(gamma > metric_rtol * max(gamma, 0))
    if (length(keep) == 0L) stop("R is numerically zero; no components can be extracted.")
    U <- ER$vectors[, keep, drop = FALSE]
    gam_sqrt_inv <- 1 / sqrt(gamma[keep])
  }

  # Work in reduced coordinates: S = Lambda^{-1/2} (U^T C U) Lambda^{-1/2},
  # applied as row/column scaling rather than dense diagonal matmuls.
  CU <- as.matrix(Matrix::crossprod(U, C %*% U))   # U^T C U
  Sred <- gam_sqrt_inv * CU
  Sred <- sweep(Sred, 2, gam_sqrt_inv, `*`)
  Sred <- 0.5 * (Sred + t(Sred))
  ES <- eigen(Sred, symmetric = TRUE)

  # Filter positive eigenvalues
  lam_all <- pmax(ES$values, 0)
  pos <- which(lam_all > 0 & lam_all > rank_rtol^2 * max(lam_all, 0))
  if (length(pos) == 0L) stop("No positive eigenvalues found (within tol).")
  if (is.null(ncomp)) ncomp <- length(pos)
  ncomp <- min(ncomp, length(pos))
  lam <- lam_all[pos][1:ncomp]
  W  <- ES$vectors[, pos, drop = FALSE][, 1:ncomp, drop = FALSE]

  # Map back: V = U Lambda^{-1/2} W  (R-orthonormal: V^T R V = I)
  V <- as.matrix(U %*% (gam_sqrt_inv * W))

  # Explained variance: sum(lambda) equals trace(R^{-1/2} C R^{-1/2}) over range(R)
  total_var <- sum(diag(Sred))
  propv <- as.numeric(lam / ifelse(total_var > 0, total_var, 1))
  cumv  <- cumsum(propv)

  list(
    v      = Matrix::Matrix(V, sparse = FALSE),  # loadings (p x k), R-orthonormal
    d      = sqrt(lam),                          # singular values (sqrt of variances)
    lambda = lam,                                # variances under the R-metric
    k      = ncomp,                              # number of components returned
    propv  = propv,
    cumv   = cumv,
    R_rank = length(keep),
    method = "geigen"
  )
}
