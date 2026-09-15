#' @title Prepare and validate constraint matrices
#' @description Coerces `A`/`M` (NULL, weight vector, diagonal, dense or
#' sparse matrix) to `Matrix` objects and validates them: finite entries,
#' symmetry within roundoff (see [symmetrize_or_stop()]), and positive
#' semi-definiteness within `tol` relative to the scale of the matrix. The
#' requested `remedy` is applied to a metric that fails the PSD check (or has
#' negative eigenvalues within tolerance when explicit clipping is requested),
#' and every repair emits a warning of class `genpca_metric_repaired` carrying
#' the [repair_metric()] report; valid PSD metrics, singular ones included,
#' pass through under every remedy. Explicit clipping removes even negative
#' eigenvalues within the validation tolerance. Asymmetric input is an error
#' under every remedy.
#' @param X data matrix (only its dimensions are used)
#' @param A,M column/row constraints
#' @param tol relative PSD tolerance (default `sqrt(.Machine$double.eps)`)
#' @param remedy what to do with an indefinite metric
#' @param verbose emit a message when a metric is replaced by the identity
#' @return list with elements `A` and `M`
#' @keywords internal
#' @importFrom assertthat assert_that
prep_constraints <- function(X, A, M, tol = .metric_rtol_default(),
                             remedy = c("error", "ridge", "clip", "identity"),
                             verbose = FALSE) {
  remedy <- match.arg(remedy)
  list(A = .prep_one_metric(A, ncol(X), "A", tol, remedy, verbose),
       M = .prep_one_metric(M, nrow(X), "M", tol, remedy, verbose))
}

.prep_diag_metric <- function(w, dim, name, rtol, remedy, verbose, from_vector = FALSE) {
  w <- as.numeric(w)
  if (any(!is.finite(w))) {
    stop("Diagonal elements of ", name, " must be finite", call. = FALSE)
  }
  if (remedy == "clip" && any(w < 0)) {
    B <- repair_metric(Matrix::Diagonal(dim, x = w), method = "clip", rtol = rtol, name = name)
    .warn_metric_repaired(attr(B, "repair_report"))
    attr(B, "repair_report") <- NULL
    return(B)
  }
  s <- max(abs(w), 0)
  if (s > 0 && any(w < -rtol * s)) {
    if (remedy == "error") {
      if (from_vector) {
        stop("Diagonal elements of ", name, " (from vector) must be non-negative (smallest = ",
             signif(min(w), 3), ")", call. = FALSE)
      }
      stop("Matrix ", name, " must be positive semi-definite (negative diagonal element ",
           signif(min(w), 3), ")", call. = FALSE)
    }
    if (verbose && remedy == "identity") message("Matrix ", name, " is not SPD, replacing with identity matrix")
    B <- repair_metric(Matrix::Diagonal(dim, x = w), method = remedy, rtol = rtol, name = name)
    .warn_metric_repaired(attr(B, "repair_report"))
    attr(B, "repair_report") <- NULL
    return(B)
  } else if (any(w < 0)) {
    message("Setting ", sum(w < 0), " tiny negative weight(s) in ", name,
            " to zero (all within ", signif(rtol, 3), " * max weight).")
    w[w < 0] <- 0
  }
  Matrix::Diagonal(n = dim, x = w)
}

.prep_one_metric <- function(W, dim, name, rtol, remedy, verbose) {
  if (is.null(W)) return(Matrix::Diagonal(dim))
  other <- if (name == "A") "ncol(X)" else "nrow(X)"
  if (is.numeric(W) && is.null(dim(W))) {
    assert_that(length(W) == dim, msg = paste0("Length of vector ", name, " must equal ", other))
    return(.prep_diag_metric(W, dim, name, rtol, remedy, verbose, from_vector = TRUE))
  }
  if (!methods::is(W, "Matrix")) W <- Matrix::Matrix(W, sparse = FALSE)
  assert_that(nrow(W) == dim, msg = paste("nrow(", name, ") != ", other, " -- ", nrow(W), " != ", dim))
  assert_that(ncol(W) == dim, msg = paste("ncol(", name, ") != ", other, " -- ", ncol(W), " != ", dim))
  if (Matrix::isDiagonal(W)) {
    return(.prep_diag_metric(Matrix::diag(W), dim, name, rtol, remedy, verbose))
  }
  xvals <- if (methods::is(W, "sparseMatrix")) W@x else as.numeric(as.matrix(W))
  if (length(xvals) && !all(is.finite(xvals))) {
    stop("Matrix ", name, " contains non-finite values", call. = FALSE)
  }
  W <- symmetrize_or_stop(W, name = name)
  if (remedy == "clip" || !is_psd(W, rtol = rtol)) {
    if (remedy == "error") {
      stop("Matrix ", name, " must be positive semi-definite", call. = FALSE)
    }
    if (verbose && remedy == "identity") message("Matrix ", name, " is not SPD, replacing with identity matrix")
    W <- repair_metric(W, method = remedy, rtol = rtol, name = name)
    if (isTRUE(attr(W, "repair_report")$changed)) .warn_metric_repaired(attr(W, "repair_report"))
    attr(W, "repair_report") <- NULL
  }
  .standardize_metric(W)
}

# Diagonal and sparse stay as they are; dense symmetric classes become
# general dense (dgeMatrix), which every downstream product accepts.
.standardize_metric <- function(W) {
  if (methods::is(W, "ddiMatrix") || methods::is(W, "sparseMatrix")) return(W)
  if (methods::is(W, "dsyMatrix") || methods::is(W, "dpoMatrix")) return(as_dge(W))
  W
}

# TRUE when a multivarious pre-processor is a single pass() step, so sparse X
# can skip the dense fit_transform path.
is_pass_preproc <- function(preproc) {
  steps <- NULL
  if (inherits(preproc, "prepper")) {
    steps <- preproc$steps
  } else if (inherits(preproc, "pre_processor") && !is.null(preproc$preproc$steps)) {
    steps <- preproc$preproc$steps
  }

  length(steps) == 1L && inherits(steps[[1L]], "pass")
}

#' Generalised Principal Components Analysis (GPCA)
#'
#' Implements the Generalised Least-Squares Matrix Decomposition of
#' Allen, Grosenick & Taylor (2014) for data observed in a **row**
#' inner-product space M and a **column** inner-product space A.
#' Setting M = I_n, A = I_p recovers ordinary PCA.
#'
#' @section Method:
#' We compute the rank-ncomp factors UDVT that minimise
#' \deqn{ \|X - UDV^\top\|_{M,A}^2
#'       = \mathrm{tr}\!\bigl(M\, (X-UDV^\top)\,A\,(X-UDV^\top)^\top\bigr) }
#' subject to UT M U = I, VT AV = I. (Allen et al., 2014).
#' Five methods are available via the `method` argument:
#' \itemize{
#'  \item{\code{"eigen"} (Default): Uses a one-shot eigen decomposition strategy based on \code{gmdLA}. It explicitly forms and decomposes a \eqn{p \times p} or \eqn{n \times n} matrix (depending on \code{n} vs \code{p}).}
#'  \item{\code{"auto"}: Chooses among \code{"eigen"}, \code{"spectra"}, and \code{"randomized"} using heuristics on shape, rank ratio (\code{ncomp / min(n,p)}), and constraint structure.}
#'  \item{\code{"spectra"}: Computes the top-k singular triplets of the metric-whitened data \eqn{F_M' X F_A} (with \eqn{M = F_M F_M'}, \eqn{A = F_A F_A'}) as an implicit operator via the \pkg{eigencore} package, without forming the large intermediate matrix. Generally faster and uses less memory for large \code{n} or \code{p} when few components are requested.}
#'  \item{\code{"randomized"}: Uses a randomized block range finder and small projected eigendecomposition. This is an approximate low-pass method that is often much faster for wide dense matrices with sparse metrics when only top components are needed.}
#'  \item{\code{"deflation"}: Uses an iterative power/deflation algorithm. Can be slower but potentially uses less memory than \code{"eigen"} for very large dense problems where \code{ncomp} is small.}
#' }
#'
#' @section Backend Guidance:
#' The default is \code{method = "eigen"}; \code{"auto"} is opt-in, not the
#' default.
#' \itemize{
#'   \item Use \code{"eigen"} (the default) when you need a stable reference solution on small/medium problems.
#'   \item Use \code{"auto"} to let a heuristic pick among \code{"eigen"}, \code{"spectra"}, and \code{"randomized"} based on problem shape and constraint structure.
#'   \item Use \code{"spectra"} for larger matrix-free iterative solves where memory pressure is a concern.
#'   \item Use \code{"randomized"} for wide low-rank settings (\code{p >> n}) with sparse metrics when throughput matters most.
#'   \item Use \code{"deflation"} when you only need a few components and can tolerate iterative convergence behavior.
#' }
#'
#' For pre-computed covariance matrices C = X'MX, see \code{\link{genpca_cov}} which
#' performs GPCA directly on C with column constraint R (equivalent to A).
#'
#' @param X   Numeric matrix n x p.
#' @param A   Column constraint: vector (implies diagonal), dense matrix, or sparse
#'            symmetric p x p PSD matrix. If `NULL`, defaults to identity.
#' @param M   Row constraint: vector (implies diagonal), dense matrix, or sparse
#'            symmetric n x n PSD matrix. If `NULL`, defaults to identity.
#' @param ncomp Number of components to extract. Defaults to `min(dim(X))`. Must be positive.
#' @param method Character string specifying the computation method. One of \code{"eigen"} (default, uses \code{gmdLA}), \code{"auto"} (heuristic choice among \code{"eigen"}, \code{"spectra"}, and \code{"randomized"}), \code{"spectra"} (iterative partial SVD of the metric-whitened data via \pkg{eigencore}, \code{gmd_spectra}), \code{"randomized"} (approximate randomized block solver \code{gmd_randomized}), or \code{"deflation"} (uses \code{gmd_deflationR} or \code{gmd_deflation_cpp}).
#' @param constraints_remedy Character string specifying what to do with a
#'        supplied `A` or `M` that is not positive semi-definite (within a
#'        relative tolerance of `sqrt(.Machine$double.eps)`). Default
#'        `"error"`: reject the input. The alternatives repair it and emit a
#'        warning of class `genpca_metric_repaired` whose `report` field
#'        (see \code{\link{repair_metric}}) records the minimum eigenvalue
#'        before and after, the shift applied, the rank and the condition
#'        number: \code{"ridge"} (Gershgorin diagonal shift: add the smallest
#'        diagonal loading that restores positive definiteness, falling back
#'        to `Matrix::nearPD()` for small dense matrices), \code{"clip"}
#'        (spectral clip to the PSD cone by zeroing negative eigenvalues;
#'        this densifies the matrix and refuses sparse input larger than
#'        2000 rows/cols, where \code{"ridge"} should be used instead), or
#'        \code{"identity"} (replace the matrix with the identity). An
#'        asymmetric metric is an error under every setting. Singular PSD
#'        metrics are valid input and are never repaired.
#' @param preproc Pre-processing transformer object from the **multivarious** package
#'                (default `multivarious::pass()`). Use `multivarious::center()` for centered GPCA.
#'                See `?multivarious::prep` for options.
#' @param threshold Convergence tolerance for the \code{"deflation"} method's
#'        inner loop. Default `1e-6`. Cutoffs are relative to the scale of
#'        the problem (the norm/singular-value floors scale with
#'        \eqn{\sqrt{\mathrm{tr}(X'MXA)}}), so results are invariant to
#'        rescaling `X`. The convergence check is on a squared step
#'        difference, so the resulting singular-vector accuracy scales like
#'        \eqn{\sqrt{\code{threshold}}}, not \code{threshold} itself.
#' @param maxit_deflation Maximum iterations per component for the
#'        \code{"deflation"} method. Default `500`.
#' @param use_cpp Logical. If `TRUE` (default) and package was compiled with C++ support,
#'                use faster C++ implementation for \code{method = "deflation"}. Fallback to R otherwise.
#'                (Ignored for \code{method = "eigen"} and \code{method = "spectra"}).
#' @param maxeig For \code{method = "eigen"} and \code{method = "spectra"}: a
#'               positive definite general metric is factored exactly by
#'               Cholesky at any size, but a singular general metric (e.g. a
#'               graph Laplacian) needs a dense eigendecomposition of the
#'               metric, which is refused when the metric has more than
#'               \code{maxeig} rows. The error names the alternatives
#'               (\code{method = "deflation"}, which only multiplies by the
#'               metric, or raising \code{maxeig}); \code{method = "auto"}
#'               routes such cases to deflation. Results are never
#'               approximated. Default `5000`.
#' @param warn_approx Deprecated and ignored: \code{method = "eigen"} no
#'        longer approximates anything.
#' @param maxit_spectra Retained for compatibility and currently unused: the
#'        \pkg{eigencore} partial SVD used by \code{method = "spectra"} is
#'        controlled by \code{tol_spectra} alone.
#' @param tol_spectra Convergence tolerance of the iterative solver when
#'        \code{method = "spectra"}. Default `1e-9`. This governs iteration
#'        only; rank decisions use \code{rank_rtol}.
#' @param rank_rtol Relative cutoff for component acceptance, on the scale of
#'        the singular values: component \eqn{j} is dropped when
#'        `d_j <= rank_rtol * d_1`. Applied by every method
#'        (for the eigen paths on \eqn{d_j^2}), so the number of components
#'        returned does not change when `X` is rescaled. Default `1e-6`.
#'        Metric validation uses a separate relative tolerance,
#'        `sqrt(.Machine$double.eps)`, for positive semi-definiteness and
#'        null-space detection for general metric eigendecompositions.
#'        Every strictly positive diagonal weight is retained in both the
#'        forward and inverse factors.
#' @param oversample Oversampling for \code{method = "randomized"} (sketch size = \code{ncomp + oversample}). Default `20`.
#' @param n_power Number of power iterations for \code{method = "randomized"}. Default `1`.
#' @param n_polish Number of optional block-polish iterations for \code{method = "randomized"}. Default `0`.
#' @param jitter_metric Relative Gram jitter for the candidate Cholesky
#'        preconditioner in \code{method = "randomized"}. The basis is checked
#'        in the original metric; a failed check uses rank-revealing
#'        orthonormalization instead. Default `1e-10`.
#' @param seed_randomized Optional seed for \code{method = "randomized"}.
#'        Default `1234`. This fully determines the randomized backend's
#'        random stream: the C++ kernel seeds its own generator from this
#'        value rather than from R's `set.seed()`/`.Random.seed`, and calling
#'        `genpca()` with `method = "randomized"` does not alter the
#'        caller's `.Random.seed`. To reproduce a randomized fit, fix
#'        `seed_randomized`, not the R seed.
#' @param tol_polish_randomized Relative tolerance used for early stopping of polish iterations in \code{method = "randomized"}. Set `0` to disable early stop. Default `1e-4`.
#' @param verbose Logical. If `TRUE`, print progress messages. Default `FALSE`.
#'
#' @return An object of class `c("genpca", "bi_projector")` inheriting from `multivarious::bi_projector`,
#'   with slots including:
#'   \describe{
#'     \item{u,v}{Left/right singular vectors scaled by the constraint metrics
#'                (MU, AV). These correspond to components in the original space's geometry.
#'                Use `components(fit)`.}
#'     \item{ou,ov}{Orthonormal singular vectors in the constraint metric
#'                  (U, V such that UT M U = I, VT AV = I). These are the core mathematical factors.}
#'     \item{sdev}{Generalised singular values d_k. Note these are singular
#'                 values of the metric-whitened data matrix, not standard
#'                 deviations: with identity metrics and centering,
#'                 `sdev = prcomp(X)$sdev * sqrt(nrow(X) - 1)`.}
#'     \item{s}{Scores: the generalised principal components
#'              `z_k = X A ov_k = ou_k d_k` (Allen et al. 2014, Section 2.4).
#'              Identical to `project(fit, X)` on the training data.
#'              Use `scores(fit)`.}
#'     \item{preproc}{The `multivarious` pre-processing object used.}
#'     \item{A, M}{The constraint matrices used (potentially after coercion to sparse format).}
#'     \item{propv}{Proportion of generalized variance explained by each component.}
#'     \item{cumv}{Cumulative proportion of generalized variance explained.}
#'   }
#'
#' @references
#' Allen, G. I., Grosenick, L., & Taylor, J. (2014).
#' *A Generalized Least-Squares Matrix Decomposition.*
#' Journal of the American Statistical Association, 109(505), 145-159.
#' arXiv:1102.3074.
#'
#' @seealso \code{\link{genpca_cov}} for GPCA on pre-computed covariance matrices,
#'   \code{\link{truncate.genpca}}, \code{\link{reconstruct.genpca}},
#'   `multivarious::bi_projector`, `multivarious::project`, `multivarious::scores`,
#'   `multivarious::components`, `multivarious::reconstruct`.
#'
#' @examples
#' if (requireNamespace("multivarious", quietly = TRUE)) {
#'   set.seed(123)
#'   X <- matrix(stats::rnorm(200 * 100), 200, 100)
#'   rownames(X) <- paste0("R", 1:200)
#'   colnames(X) <- paste0("C", 1:100)
#'
#'   # Standard PCA (A=I, M=I, centered) - using default method="eigen"
#'   gpca_std_eigen <- genpca(X, ncomp = 5, preproc = multivarious::center(), verbose = FALSE)
#'
#'   # Standard PCA using Spectra method (requires C++ build)
#'   # gpca_std_spectra <- try(genpca(X, ncomp = 5,
#'   #                              preproc = multivarious::center(),
#'   #                              method = "spectra", verbose = TRUE))
#'   # if (!inherits(gpca_std_spectra, "try-error")) {
#'   #    print(head(gpca_std_spectra$sdev))
#'   # }
#'
#'   # Compare singular values with prcomp
#'   pr_std <- stats::prcomp(X, center = TRUE, scale. = FALSE)
#'   print("Eigen Method Sdev:")
#'   print(head(gpca_std_eigen$sdev))
#'   print("prcomp Sdev:")
#'   print(head(pr_std$sdev))
#'   print(paste("Total Var Explained (Eigen):",
#'               round(sum(gpca_std_eigen$propv) * 100), "%"))
#'
#'   # Weighted column PCA (diagonal A, no centering)
#'   col_weights <- stats::runif(100, 0.5, 1.5)
#'   gpca_weighted <- genpca(X, A = col_weights, ncomp = 3,
#'                           preproc = multivarious::pass(), verbose = FALSE)
#'   print("Weighted GPCA Sdev:")
#'   print(gpca_weighted$sdev)
#'   print(head(components(gpca_weighted)))
#' }
#' @useDynLib genpca, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom multivarious bi_projector fit_transform pass scores sdev components reconstruct inverse_transform ncomp
#' @importFrom Matrix Matrix isSymmetric isDiagonal diag t forceSymmetric Diagonal crossprod tcrossprod
#' @importFrom assertthat assert_that
#' @importFrom methods as is
#' @importFrom stats rnorm runif
#' @export
genpca <- function(X, A = NULL, M = NULL, ncomp = NULL,
                   method = c("eigen", "auto", "spectra", "randomized", "deflation"),
                   constraints_remedy = c("error", "ridge", "clip", "identity"),
                   preproc = multivarious::pass(), # Default to pass() for safety
                   threshold = 1e-6, # For deflation
                   maxit_deflation = 500L, # For deflation
                   use_cpp = TRUE, # For deflation
                   maxeig = 5000, # For method="eigen": bound on dense eigendecomposition of a singular metric
                   warn_approx = TRUE, # Deprecated: no approximation is made any more
                   maxit_spectra = 1000, # For method="spectra"
                   tol_spectra = 1e-9,   # For method="spectra"
                   rank_rtol = 1e-6,     # Relative singular-value cutoff (all methods)
                   oversample = 20L,     # For method="randomized"
                   n_power = 1L,         # For method="randomized"
                   n_polish = 0L,        # For method="randomized"
                   jitter_metric = 1e-10, # For method="randomized"
                   seed_randomized = 1234L, # For method="randomized"
                   tol_polish_randomized = 1e-4, # For method="randomized"
                   verbose = FALSE) {

  method <- match.arg(method)
  constraints_remedy <- match.arg(constraints_remedy)
  if (!missing(warn_approx)) {
    warning("`warn_approx` is deprecated and ignored: method = 'eigen' no longer approximates.", call. = FALSE)
  }

  if (is.null(ncomp)) {
      ncomp <- min(dim(X))
  }
  assert_that(length(ncomp) == 1 && ncomp == floor(ncomp) && ncomp > 0,
              msg = "ncomp must be a single positive integer.")
  ncomp <- min(min(dim(X)), ncomp) # Cannot exceed dimensions
  assert_that(length(maxit_deflation) == 1 &&
                maxit_deflation == floor(maxit_deflation) &&
                maxit_deflation > 0,
              msg = "maxit_deflation must be a single positive integer.")
  maxit_deflation <- as.integer(maxit_deflation)
  assert_that(length(oversample) == 1 &&
                oversample == floor(oversample) &&
                oversample >= 0,
              msg = "oversample must be a single non-negative integer.")
  oversample <- as.integer(oversample)
  assert_that(length(n_power) == 1 &&
                n_power == floor(n_power) &&
                n_power >= 0,
              msg = "n_power must be a single non-negative integer.")
  n_power <- as.integer(n_power)
  assert_that(length(n_polish) == 1 &&
                n_polish == floor(n_polish) &&
                n_polish >= 0,
              msg = "n_polish must be a single non-negative integer.")
  n_polish <- as.integer(n_polish)
  if (!is.null(seed_randomized)) {
    assert_that(length(seed_randomized) == 1 &&
                  seed_randomized == floor(seed_randomized),
                msg = "seed_randomized must be NULL or a single integer.")
    seed_randomized <- as.integer(seed_randomized)
  }
  assert_that(length(tol_polish_randomized) == 1 &&
                is.finite(tol_polish_randomized) &&
                tol_polish_randomized >= 0,
              msg = "tol_polish_randomized must be a single non-negative number.")
  assert_that(is.numeric(rank_rtol) && length(rank_rtol) == 1 &&
                is.finite(rank_rtol) && rank_rtol >= 0,
              msg = "rank_rtol must be a single non-negative number.")

  # Prepare and validate constraints M and A
  if (verbose) message("Preparing constraints...")
  pcon <- prep_constraints(X, A, M, remedy = constraints_remedy, verbose = verbose)
  A <- pcon$A
  M <- pcon$M

  # Prepare and apply pre-processing using fit_transform API
  if (verbose) message("Applying pre-processing...")
  if (methods::is(X, "sparseMatrix") && is_pass_preproc(preproc)) {
    ft <- multivarious::fit_transform(preproc, matrix(0, nrow = 1L, ncol = ncol(X)))
    procres <- ft$preproc
    Xp <- X
  } else {
    # multivarious pre-processors require a base matrix. Any non-pass
    # transform (e.g. centering) densifies a sparse X anyway, so convert.
    if (methods::is(X, "sparseMatrix")) {
      if (verbose) message("Densifying sparse X for pre-processing (preproc is not pass()).")
      X <- as.matrix(X)
    }
    ft <- multivarious::fit_transform(preproc, X)
    procres <- ft$preproc
    Xp <- ft$transformed
  }

  n <- nrow(Xp)
  p <- ncol(Xp)

  # Check if C++ code is available (specific function name depends on package build)
  # Placeholder check - replace with actual check if package uses compiled code
  cpp_deflation_available <- exists("gmd_deflation_cpp", mode = "function") # Example check
  cpp_spectra_available <- exists("gmd_spectra", mode = "function")
  cpp_randomized_available <- exists("gmd_randomized_cpp_dn", mode = "function")

  selected_method <- method
  if (method == "auto") {
      min_dim <- min(n, p)
      max_dim <- max(n, p)
      k_ratio <- ncomp / min_dim
      diagonal_constraints <- Matrix::isDiagonal(A) && Matrix::isDiagonal(M)
      sparse_constraints <- methods::is(A, "sparseMatrix") || methods::is(M, "sparseMatrix")
      wide_problem <- p >= (2L * n)
      large_problem <- (min_dim >= 120L) && (max_dim >= 1500L)
      use_randomized <- cpp_randomized_available &&
        sparse_constraints &&
        !diagonal_constraints &&
        wide_problem &&
        large_problem &&
        k_ratio <= 0.15
      use_spectra <- cpp_spectra_available && (
        (diagonal_constraints && min_dim >= 200L && k_ratio <= 0.35) ||
        (!diagonal_constraints && min_dim >= 1000L && k_ratio <= 0.10) ||
        (sparse_constraints && min_dim >= 100L && k_ratio <= 0.25)
      )
      # A singular general metric on the small side needs a dense
      # eigendecomposition under both "eigen" and "spectra"; above maxeig only
      # deflation can proceed (it multiplies by the metric, never factors it).
      small_metric <- if (p <= n) A else M
      needs_dense_eigen <- !diagonal_constraints && min_dim > maxeig &&
        !Matrix::isDiagonal(small_metric) && !is_pd(small_metric)
      selected_method <- if (needs_dense_eigen) {
        "deflation"
      } else if (use_randomized) {
        "randomized"
      } else if (use_spectra) {
        "spectra"
      } else {
        "eigen"
      }
      if (verbose) {
        message(
          "Auto method selected '", selected_method,
          "' (min_dim=", min_dim,
          ", max_dim=", max_dim,
          ", k_ratio=", sprintf("%.3f", k_ratio),
          ", diagonal_constraints=", diagonal_constraints,
          ", sparse_constraints=", sparse_constraints,
          ", wide_problem=", wide_problem, ")."
        )
      }
  }

  if (selected_method == "deflation" && use_cpp && !cpp_deflation_available) {
      if (verbose) message("use_cpp=TRUE but C++ deflation code not found. Falling back to R version.")
      use_cpp <- FALSE # Force R version if C++ not found
  }
  if (selected_method == "spectra" && !cpp_spectra_available) {
      stop("method='spectra' requires the internal solver 'gmd_spectra', which was not found.")
  }

  # --- Core Decomposition --- #
  if (selected_method == "deflation") {
    if (verbose) message(paste0("Using iterative deflation (", ifelse(use_cpp, "C++", "R"), ") to extract ", ncomp, " components..."))
    if (use_cpp) {
      # C++ deflation backend expects sparse metrics.
      M_cpp <- as_dgc(M)
      A_cpp <- as_dgc(A)
    }
    if (n < p) {
        if (use_cpp) {
          svdfit <- gmd_deflation_cpp_dispatch(Matrix::t(Xp), A_cpp, M_cpp, ncomp,
                                               thr = threshold, rank_rtol = rank_rtol,
                                               maxit = maxit_deflation,
                                               verbose = verbose)
          if (is.matrix(svdfit$d)) {
            svdfit$d <- svdfit$d[, 1] # Ensure d is vector
          } else {
            svdfit$d <- as.vector(svdfit$d)
          }
          # C++ might not return propv/cumv, need to calculate if necessary
          if (is.null(svdfit$k)) svdfit$k <- length(svdfit$d)
        } else {
          svdfit <- gmd_deflationR(Matrix::t(Xp), A, M, ncomp,
                                   thr = threshold, rank_rtol = rank_rtol,
                                   maxit = maxit_deflation,
                                   verbose = verbose)
        }
        # Swap u and v back
        svdfit <- list(u = svdfit$v, v = svdfit$u, d = svdfit$d, k = svdfit$k, cumv = svdfit$cumv, propv = svdfit$propv)
    } else {
        if (use_cpp) {
          svdfit <- gmd_deflation_cpp_dispatch(Xp, M_cpp, A_cpp, ncomp,
                                               thr = threshold, rank_rtol = rank_rtol,
                                               maxit = maxit_deflation,
                                               verbose = verbose)
          if (is.matrix(svdfit$d)) {
            svdfit$d <- svdfit$d[, 1]
          } else {
            svdfit$d <- as.vector(svdfit$d)
          }
          if (is.null(svdfit$k)) svdfit$k <- length(svdfit$d)
        } else {
          svdfit <- gmd_deflationR(Xp, M, A, ncomp,
                                   thr = threshold, rank_rtol = rank_rtol,
                                   maxit = maxit_deflation,
                                   verbose = verbose)
        }
    }
    # Deflation methods should return propv/cumv, but double check
    if (is.null(svdfit$propv) || is.null(svdfit$cumv)) {
       if (verbose) message(" Calculating variance explained for deflation method...")
       total_variance <- sum((M %*% Xp) * (Xp %*% A)) # tr(Xp' M Xp A) without forming p x p
       if (is.finite(total_variance) && total_variance > 0) {
          svdfit$propv <- svdfit$d^2 / total_variance
          svdfit$cumv <- cumsum(svdfit$propv)
       } else {
          svdfit$propv <- rep(0, svdfit$k)
          svdfit$cumv <- rep(0, svdfit$k)
       }
    }

  } else if (selected_method == "eigen") { # One-shot eigen-decomposition approach
    if (verbose) message(paste0("Using one-shot eigen decomposition (gmdLA) to extract ", ncomp, " components..."))
    if (n < p) {
        if (verbose) message(" (n < p, using dual formulation)")
        # Dual formulation works on the n x n problem directly: eigendecompose
        # M^{1/2} (X A X') M^{1/2} instead of the p x p primal target.
        svdfit <- gmdLA(Xp, M, A, k = ncomp, n_orig = n, p_orig = p,
                        maxeig = maxeig, rank_rtol = rank_rtol, use_dual = TRUE,
                        warn_approx = warn_approx, verbose = verbose)
    } else {
        svdfit <- gmdLA(Xp, M, A, k = ncomp, n_orig = n, p_orig = p,
                        maxeig = maxeig, rank_rtol = rank_rtol, use_dual = FALSE,
                        warn_approx = warn_approx, verbose = verbose)
    }

  } else if (selected_method == "spectra") { # Matrix-free C++/Spectra approach
      if (verbose) message(paste0("Using the iterative whitened-operator SVD (eigencore) to extract ", ncomp, " components..."))
      # Ensure Xp is dense matrix for the C++ function
      Xp_dense <- as.matrix(Xp)
      if (any(!is.finite(Xp_dense))) stop("Input matrix X (after preproc) contains non-finite values.")

      # Call the C++ function
      spectra_res <- tryCatch(gmd_spectra(Xp_dense, M, A, k = ncomp, tol = tol_spectra,
                                          maxit = maxit_spectra, rank_rtol = rank_rtol,
                                          dense_maxn = maxeig),
                              error = function(e) {stop("Call to gmd_spectra failed: ", e$message)})

      # Calculate variance explained
      if (verbose) message(" Calculating variance explained for Spectra method...")
      total_variance <- sum((M %*% Xp) * (Xp %*% A)) # tr(Xp' M Xp A) without forming p x p
      if (!is.finite(total_variance) || total_variance <= 0) {
          propv <- rep(0, spectra_res$k)
          warning("Total generalized variance is near zero.")
      } else {
          propv <- (spectra_res$d^2) / total_variance
      }
      cumv <- cumsum(propv)

      # Map results to the svdfit structure
      # IMPORTANT: spectra_res returns scores and components, not orthonormal eigenvectors!
      # We need to back-calculate the orthonormal eigenvectors for consistency
      # scores = M * ou * D, so ou = M^{-1} * scores / D
      # components = A * ov, so ov = A^{-1} * components

      # For now, set a flag to handle this differently later
      svdfit <- list(d = spectra_res$d,
                     u = spectra_res$u, # These are actually scores/D, not ou!
                     v = spectra_res$v, # These are actually components, not ov!
                     ou = spectra_res$ou,
                     ov = spectra_res$ov,
                     k = spectra_res$k,
                     propv = propv,
                     cumv = cumv,
                     is_spectra = TRUE) # Flag to indicate special handling needed
  } else if (selected_method == "randomized") { # Randomized block range finder
      if (verbose) {
        message(
          "Using randomized block solver to extract ", ncomp,
          " components (oversample=", oversample,
          ", n_power=", n_power,
          ", n_polish=", n_polish,
          ", tol_polish=", tol_polish_randomized, ")."
        )
      }
      Xp_dense <- as.matrix(Xp)
      if (any(!is.finite(Xp_dense))) stop("Input matrix X (after preproc) contains non-finite values.")

      rand_res <- tryCatch(
        gmd_randomized(
          X = Xp_dense,
          Q = M,
          R = A,
          k = ncomp,
          oversample = oversample,
          n_power = n_power,
          n_polish = n_polish,
          jitter = jitter_metric,
          tol = rank_rtol,
          polish_tol = tol_polish_randomized,
          seed = seed_randomized
        ),
        error = function(e) stop("Call to gmd_randomized failed: ", e$message)
      )

      if (verbose) message(" Calculating variance explained for randomized method...")
      total_variance <- sum((M %*% Xp) * (Xp %*% A)) # tr(Xp' M Xp A) without forming p x p
      if (!is.finite(total_variance) || total_variance <= 0) {
        propv <- rep(0, rand_res$k)
        warning("Total generalized variance is near zero.")
      } else {
        propv <- (rand_res$d^2) / total_variance
      }
      cumv <- cumsum(propv)

      svdfit <- list(
        d = rand_res$d,
        u = rand_res$u, # Q-orthonormal vectors
        v = rand_res$v, # R-orthonormal vectors
        k = rand_res$k,
        propv = propv,
        cumv = cumv
      )
  } else {
      stop("Internal error: Unknown method specified.") # Should not happen due to match.arg
  }

  # Enforce the public component cutoff at the common boundary, including
  # every vector and variance summary returned by a backend.
  keep <- .keep_components(svdfit$d, rank_rtol)
  for (nm in c("u", "v", "ou", "ov")) {
    if (!is.null(svdfit[[nm]])) svdfit[[nm]] <- svdfit[[nm]][, keep, drop = FALSE]
  }
  svdfit$d <- svdfit$d[keep]
  svdfit$propv <- svdfit$propv[keep]
  svdfit$cumv <- cumsum(svdfit$propv)
  svdfit$k <- length(svdfit$d)

  # Check how many components were actually found
  k_found <- svdfit$k # gmdLA and gmd_deflationR/cpp should return 'k'
  if (is.null(k_found)) k_found <- length(svdfit$d) # Fallback if k not returned

  if (k_found < ncomp) {
      warning("Found only ", k_found, " valid components, less than requested ncomp=", ncomp)
      ncomp <- k_found
      # Trim results if necessary (should be done in helpers, but ensures consistency)
      if (length(svdfit$d) > ncomp) svdfit$d <- svdfit$d[1:ncomp]
      if (ncol(svdfit$u) > ncomp) svdfit$u <- svdfit$u[, 1:ncomp, drop = FALSE]
      if (ncol(svdfit$v) > ncomp) svdfit$v <- svdfit$v[, 1:ncomp, drop = FALSE]
      if (length(svdfit$propv) > ncomp) svdfit$propv <- svdfit$propv[1:ncomp]
      if (length(svdfit$cumv) > ncomp) svdfit$cumv <- svdfit$cumv[1:ncomp]
  }

  if (ncomp == 0) {
      warning("No valid components found.")
      # Return an empty but valid structure
      return(multivarious::bi_projector(v = matrix(0.0, p, 0), s = matrix(0.0, n, 0), sdev = numeric(0),
                                        preproc = procres, ov = matrix(0.0, p, 0), ou = matrix(0.0, n, 0),
                                        u = matrix(0.0, n, 0), classes = "genpca", A = A, M = M,
                                        propv = numeric(0), cumv = numeric(0)))
  }

  # --- Construct bi_projector object --- #
  if (verbose) message("Constructing final object...")

  # Scores: z = X A ov = ou D (Allen et al. 2014, Section 2.4: z_k = X R v_k).
  # This matches project(): X (A ov) = ou D since ov'A ov = I.

  if (!is.null(svdfit$is_spectra) && svdfit$is_spectra) {
    # From gmd_fast_cpp:
    # - svdfit$ou / svdfit$ov are the metric-orthonormal factors
    # - svdfit$u contains M-weighted scores (M ou D), svdfit$v contains A ov
    if (is.null(svdfit$ou) || is.null(svdfit$ov)) {
      # Older compiled backends returned only M ou D / A ov, from which ou/ov
      # cannot be recovered without metric solves; refuse rather than return
      # silently wrong factors.
      stop("Internal error: spectra backend did not return metric-orthonormal ",
           "factors (ou/ov). Rebuild the package.")
    }
    ou <- as.matrix(svdfit$ou)
    ov <- as.matrix(svdfit$ov)
  } else {
    # Standard path for eigen/deflation methods
    ou <- as.matrix(svdfit$u)  # Ensure regular matrix
    ov <- as.matrix(svdfit$v)  # Ensure regular matrix
  }
  scores_mat <- sweep(ou, 2, svdfit$d, `*`)

  # Assign row/col names if available from original X
  # Get original indices if preproc modified them
  row_indices <- if (!is.null(attr(Xp, "row_indices"))) attr(Xp, "row_indices") else 1:nrow(X)
  col_indices <- if (!is.null(attr(Xp, "col_indices"))) attr(Xp, "col_indices") else 1:ncol(X)

  if (!is.null(rownames(X))) {
      rownames(scores_mat) <- rownames(X)[row_indices]
  } else {
      rownames(scores_mat) <- paste0("Obs", 1:n)
  }
  colnames(scores_mat) <- paste0("PC", 1:ncomp)

  # Loadings (components): v = A V (where V is ov)
  if (!is.null(svdfit$is_spectra) && svdfit$is_spectra) {
    # spectra_res$v are already R-weighted components
    loadings_mat <- as.matrix(svdfit$v)  # Ensure regular matrix
  } else {
    loadings_mat <- as.matrix(A %*% svdfit$v)  # Ensure regular matrix
  }
  if (!is.null(colnames(X))) {
      rownames(loadings_mat) <- colnames(X)[col_indices]
  } else {
      rownames(loadings_mat) <- paste0("Var", 1:p)
  }
  colnames(loadings_mat) <- paste0("PC", 1:ncomp)

  # Create the S3 object using the multivarious constructor
  M_ou <- as.matrix(M %*% ou)

  ret <- multivarious::bi_projector(
    v = loadings_mat,     # Loadings = A %*% ov
    s = scores_mat,       # Scores = ou %*% D = X %*% A %*% ov (paper's z_k)
    sdev = svdfit$d,      # Singular values
    preproc = procres,    # Preprocessing object
    ov = ov,              # Orthonormal V in A metric (ov' A ov = I)
    ou = ou,              # Orthonormal U in M metric (ou' M ou = I)
    u = M_ou,             # Metric-weighted factor M %*% ou
    classes = "genpca",   # Specific class first
    A = A,                # Store constraint matrices
    M = M,
    propv = svdfit$propv, # Proportion of variance
    cumv = svdfit$cumv   # Cumulative variance
  )
  ret$method <- selected_method
  ret$requested_method <- method

  # bi_projector should handle adding "bi_projector", "projector" classes.

  if (verbose) message("GPCA finished.")
  return(ret)
}



#' Experimental penalized-ML estimation of GPCA metrics
#'
#' Alternates between GPCA factor estimation and penalized maximum-likelihood
#' updates of the row/column metric matrices (M, A) under a Gaussian
#' matrix-normal error model with a low-rank mean. Each iteration performs
#' three exact block minimizations of a single penalized objective:
#' \enumerate{
#'   \item Fit GPCA with current \code{A}, \code{M}: the GMD theorem makes
#'         this the best rank-\code{ncomp} fit in the (M, A) norm, so it
#'         exactly minimizes the residual term.
#'   \item Update \eqn{\Sigma_r = E A E^T / p + \lambda I} (the exact block
#'         minimizer under the ridge penalty), set \code{M = solve(Sigma_r)}.
#'   \item Update \eqn{\Sigma_c = E^T M E / n + \lambda I} using the
#'         \emph{updated} \code{M} (sequential flip-flop, Dutilleul 1999),
#'         set \code{A = solve(Sigma_c)}.
#' }
#'
#' @details
#' The objective is the matrix-normal log-likelihood with a low-rank mean and
#' an inverse-Wishart-style ridge penalty
#' \eqn{\lambda\,(p\,\mathrm{tr}\,\Sigma_r^{-1} + n\,\mathrm{tr}\,\Sigma_c^{-1})}
#' (a MAP estimate). Because every block update is an exact minimizer of this
#' one objective, \code{loglik_path} is monotone non-decreasing up to
#' numerical noise. The penalty also resolves the \eqn{c\,\Sigma_r,
#' \Sigma_c/c} scale indeterminacy, so the converged metrics are the
#' penalized optimum and no rescaling is needed (\code{scale_fix = "none"},
#' the default). \code{scale_fix = "trace"} or \code{"det"} additionally
#' applies a joint reciprocal rescale at exit (row covariance normalized,
#' factor absorbed into the column covariance). The unpenalized
#' matrix-normal likelihood is invariant to that rescale, but the penalty
#' \eqn{p\lambda\,\mathrm{tr}(M) + n\lambda\,\mathrm{tr}(A)} is not, so
#' the rescaled metrics are no longer the penalized optimum; the returned
#' \code{loglik} is always evaluated at the returned metrics and
#' \code{loglik_rescale_delta} reports how far the rescale moved it. The
#' algorithm stops when the relative change in the penalized log-likelihood
#' falls below \code{tol} or \code{max_iter} is reached. Increase
#' \code{lambda} or reduce \code{ncomp} if iterations become unstable. The
#' objective is not identifiable with \code{lambda = 0}.
#'
#' @param X Numeric matrix (n x p).
#' @param ncomp Rank to extract at each GPCA step.
#' @param max_iter Maximum outer alternations (default 20).
#' @param lambda Ridge penalty weight (default 1e-3). Part of the objective
#'        (MAP interpretation), not just a numerical safeguard: it shrinks
#'        both covariances toward a multiple of the identity and pins the
#'        row/column scale split during iteration. Must be non-negative;
#'        with \code{lambda = 0} the objective loses strict convexity in the
#'        scale direction and covariances may become singular.
#' @param scale_fix Optional post-hoc reparameterization of the
#'        \code{c * Sigma_r, Sigma_c / c} split at exit. One of \code{"none"}
#'        (default: keep the penalized optimum), \code{"trace"} (row
#'        covariance scaled to mean diagonal 1) or \code{"det"} (row
#'        covariance scaled to determinant 1). Applied as a joint reciprocal
#'        rescale, so the fitted covariance \eqn{\Sigma_r \otimes \Sigma_c}
#'        and the unpenalized likelihood are unchanged, but the penalized
#'        objective generally decreases; see Details.
#' @param tol Relative tolerance on successive penalized log-likelihood change
#'        (default 1e-4) for early stopping.
#' @param method GPCA method passed to \code{genpca} (defaults to "eigen").
#' @param constraints_remedy Passed to \code{genpca}; defaults to "error".
#'        The learned metrics are inverses of positive definite matrices, so
#'        no repair fires in practice.
#' @param preproc Pre-processing transformer; defaults to \code{multivarious::pass()}.
#' @param verbose Logical; if TRUE, prints iteration diagnostics.
#' @param ... Additional arguments forwarded to \code{genpca}.
#'
#' @return A list with elements \code{fit} (a \code{genpca} fit computed with
#'         the returned metrics), \code{A}, \code{M} (learned SPD metrics),
#'         \code{loglik} (the penalized log-likelihood evaluated at the
#'         returned \code{M}, \code{A} and \code{fit}),
#'         \code{loglik_unpenalized} (the same without the \code{lambda}
#'         penalty), \code{loglik_rescale_delta} (the change in the penalty
#'         contribution caused solely by reciprocal metric rescaling; exactly
#'         zero for \code{scale_fix = "none"}), \code{loglik_refit_delta}
#'         (the remaining change from the last path value to \code{loglik},
#'         including final refitting and numerical objective reevaluation), and
#'         \code{loglik_path} (the penalized log-likelihood after each outer
#'         iteration; monotone non-decreasing up to numerical noise, since
#'         every block update exactly minimizes the shared penalized
#'         objective). Values omit additive constants and include the
#'         \code{lambda} penalty, so they are comparable across iterations
#'         and across runs with the same \code{lambda}, but not across
#'         different \code{lambda} values.
#'
#' @examples
#' if (requireNamespace("multivarious", quietly = TRUE)) {
#'   set.seed(123)
#'   X <- matrix(rnorm(40), 8, 5)
#'   res <- gpca_mle(X, ncomp = 2, max_iter = 5, lambda = 1e-3,
#'                   scale_fix = "trace", verbose = FALSE)
#'   # Learned metrics are SPD and match dimensions
#'   dim(res$A); dim(res$M)
#'   res$loglik_path
#' }
#'
#' @references
#' Dutilleul, P. (1999). *The MLE algorithm for the matrix normal distribution*.
#' Journal of Statistical Computation and Simulation, 64(2), 105-123.
#'
#' @importFrom utils tail
#' @export
gpca_mle <- function(X, ncomp = min(dim(X)),
                     max_iter = 20,
                     lambda = 1e-3,
                     scale_fix = c("none", "trace", "det"),
                     tol = 1e-4,
                     method = "eigen",
                     constraints_remedy = "error",
                     preproc = multivarious::pass(),
                     verbose = FALSE,
                     ...) {

  scale_fix <- match.arg(scale_fix)
  n <- nrow(X)
  p <- ncol(X)
  stopifnot(is.numeric(lambda), length(lambda) == 1, lambda >= 0)
  if (lambda == 0) {
    warning("gpca_mle: with lambda = 0 the objective is not identifiable in the scale ",
            "of (M, A) and the covariance updates may become singular.", call. = FALSE)
  }

  # initialise metrics as identity
  A <- Matrix::Diagonal(p)
  M <- Matrix::Diagonal(n)

  # Penalized log-likelihood (MAP objective, constants dropped):
  #   -2*ll = p*logdet(Sigma_r) + n*logdet(Sigma_c) + tr(M E A E')
  #           + p*lambda*tr(M) + n*lambda*tr(A)
  # Each block update below is the EXACT minimizer of this objective in its
  # block, so ll is monotone non-decreasing (up to numerical noise):
  #  - Xhat: GPCA/GMD gives the best rank-ncomp fit in the (M, A) norm,
  #    which minimizes the tr(M E A E') term.
  #  - Sigma_r: argmin of p*logdet(S) + tr(S^{-1}(E A E' + p*lambda*I)) is
  #    S = E A E'/p + lambda*I.
  #  - Sigma_c (with the UPDATED M -- sequential flip-flop, Dutilleul 1999):
  #    S = E' M E/n + lambda*I.
  # The lambda penalty also pins the c*Sigma_r, Sigma_c/c scale
  # indeterminacy, so no rescaling is needed inside the loop.
  pen_loglik <- function(E, M, A, logdet_r, logdet_c) {
    quad <- sum(E * as.matrix(M %*% E %*% A))
    pen <- p * lambda * sum(Matrix::diag(M)) + n * lambda * sum(Matrix::diag(A))
    -0.5 * (p * logdet_r + n * logdet_c + quad + pen)
  }

  loglik_path <- numeric(0)
  last_ll <- -Inf

  for (it in seq_len(max_iter)) {
    if (verbose) message("[gpca_mle] Iteration ", it, ": fitting GPCA...")

    fit <- genpca(X, A = A, M = M, ncomp = ncomp,
                  method = method,
                  constraints_remedy = constraints_remedy,
                  preproc = preproc,
                  verbose = verbose && it == 1, # avoid chatter each loop
                  ...)

    Xhat <- multivarious::reconstruct(fit)  # back to data scale
    E <- X - Xhat

    # Sequential (flip-flop) covariance updates: Sigma_c uses the NEW M.
    Sigma_r <- (E %*% A %*% Matrix::t(E)) / p + lambda * Matrix::Diagonal(n)
    Sigma_r <- .mle_symmetric_pd(Sigma_r)
    M <- tryCatch(Matrix::solve(Sigma_r), error = function(e) stop("Failed to invert Sigma_r: ", e$message))

    Sigma_c <- (Matrix::t(E) %*% M %*% E) / n + lambda * Matrix::Diagonal(p)
    Sigma_c <- .mle_symmetric_pd(Sigma_c)
    A <- tryCatch(Matrix::solve(Sigma_c), error = function(e) stop("Failed to invert Sigma_c: ", e$message))

    logdet_r <- as.numeric(Matrix::determinant(Sigma_r, logarithm = TRUE)$modulus)
    logdet_c <- as.numeric(Matrix::determinant(Sigma_c, logarithm = TRUE)$modulus)
    ll <- pen_loglik(E, M, A, logdet_r, logdet_c)

    loglik_path <- c(loglik_path, ll)
    if (verbose) message(sprintf("[gpca_mle]  penalized loglik = %.4f", ll))

    if (it > 1) {
      rel_change <- abs(ll - last_ll) / (abs(last_ll) + 1e-9)
      if (rel_change < tol) {
        if (verbose) message("[gpca_mle] Converged: relative loglik change < tol")
        break
      }
    }
    last_ll <- ll
  }

  # Optional post-hoc reparameterization (c*Sigma_r, Sigma_c/c): the
  # unpenalized matrix-normal likelihood is invariant to it, but the ridge
  # penalty p*lambda*tr(M) + n*lambda*tr(A) is NOT, so the returned metrics
  # are then no longer the penalized optimum. The returned `loglik` is
  # therefore always recomputed at the metrics actually returned, and the
  # size of the move is reported in `loglik_rescale_delta`.
  loglik_at_exit <- tail(loglik_path, 1)
  loglik_rescale_delta <- 0
  if (scale_fix != "none") {
    s <- if (scale_fix == "trace") {
      sum(Matrix::diag(Sigma_r)) / n
    } else { # "det"
      logdet_r <- as.numeric(Matrix::determinant(Sigma_r, logarithm = TRUE)$modulus)
      exp(logdet_r / n)
    }
    if (is.finite(s) && s > 0) {
      # The quadratic and log-determinant terms are invariant to reciprocal
      # scaling. Evaluate only the penalty change, without conflating it
      # with the final refit or inversion roundoff in objective reevaluation.
      loglik_rescale_delta <- -0.5 * lambda * (
        p * (s - 1) * sum(Matrix::diag(M)) +
        n * (1 / s - 1) * sum(Matrix::diag(A)))
      Sigma_r <- Sigma_r / s
      Sigma_c <- Sigma_c * s
      M <- M * s
      A <- A / s
    }
  }

  # Refit once so the returned fit corresponds to the returned metrics.
  fit <- genpca(X, A = A, M = M, ncomp = ncomp,
                method = method,
                constraints_remedy = constraints_remedy,
                preproc = preproc,
                verbose = FALSE,
                ...)

  # Objective at the returned (M, A, fit). Preserve any final refit and
  # numerical reevaluation discrepancy separately from the rescaling effect.
  E_final <- X - multivarious::reconstruct(fit)
  Sigma_r_final <- Matrix::solve(M)
  Sigma_c_final <- Matrix::solve(A)
  logdet_r_final <- as.numeric(Matrix::determinant(Sigma_r_final, logarithm = TRUE)$modulus)
  logdet_c_final <- as.numeric(Matrix::determinant(Sigma_c_final, logarithm = TRUE)$modulus)
  loglik_final <- pen_loglik(E_final, M, A, logdet_r_final, logdet_c_final)
  pen_final <- p * lambda * sum(Matrix::diag(M)) + n * lambda * sum(Matrix::diag(A))

  list(fit = fit,
       A = A,
       M = M,
       loglik = loglik_final,
       loglik_unpenalized = loglik_final + 0.5 * pen_final,
       loglik_rescale_delta = loglik_rescale_delta,
       loglik_refit_delta = loglik_final - loglik_at_exit - loglik_rescale_delta,
       loglik_path = loglik_path)
}



# The block covariance updates in gpca_mle() are positive definite by
# construction when lambda > 0 (Gram matrix plus lambda * I); their only
# defect is the roundoff asymmetry of the matrix products, which
# prep_constraints() would reject. Average the triangles and repair only when
# the matrix is genuinely not positive definite (lambda = 0). A tolerant
# ensure_spd() must NOT be used here: its relative PD margin fires once the
# covariance scale exceeds ~1e6 * lambda and would replace the learned metric
# by a near-identity.
.mle_symmetric_pd <- function(S) {
  S <- Matrix::forceSymmetric((S + Matrix::t(S)) / 2)
  if (is_pd(S, rtol = 0)) S else ensure_spd(S, tol = .metric_rtol_default())
}

#' @noRd
#' @importFrom Matrix Diagonal t crossprod tcrossprod diag solve isDiagonal Matrix
#' @importFrom methods as is
# gmdLA caches the factorization of the small-side metric as an attribute on
# the matrix and returns the annotated matrix so callers can reassign it
# (e.g. R <- R_fac$matrix) and reuse the factor in subsequent calls.
gmdLA <- function(X, Q, R, k = min(n_orig, p_orig), n_orig, p_orig,
                  maxeig = 5000, rank_rtol = 1e-6,
                  metric_rtol = .metric_rtol_default(), use_dual = FALSE,
                  warn_approx = TRUE, verbose = FALSE) {

  cache_attr_name <- "eigen_decomp_cache"
  # Rank decisions are relative: metric eigenvalues below metric_rtol times
  # the largest are null space; target eigenvalues below rank_rtol^2 times
  # the largest are dropped.

  # Factor the small-side metric once: A = F F' (Cholesky when positive
  # definite, eigen factor on the range otherwise; see .metric_factor). The
  # target F' (X'QX) F is similar to (X'QX) R and V = F^{-T} Z is
  # R-orthonormal, so nothing is ever truncated. `maxeig` only bounds the
  # dense eigendecomposition needed for a singular general metric.
  metric_factor_cached <- function(M, cache_attr, mat_name) {
    cached <- attr(M, cache_attr)
    if (!is.null(cached)) {
      if (verbose) message(paste(" Using cached decomposition for matrix", mat_name))
      fac <- if (!is.null(cached$apply_t)) {
        cached
      } else {
        # legacy cache entries: symmetric square root and its pseudo-inverse
        sq <- cached$sqrtm
        isq <- cached$invsqrtm
        list(kind = "legacy", ncol = ncol(sq), mat = sq,
             apply = function(V) sq %*% V,
             apply_t = function(U) sq %*% U,
             solve_t = function(Z) isq %*% Z)
      }
    } else {
      if (verbose) message(paste(" Computing factorization for matrix", mat_name))
      fac <- .metric_factor(M, metric_rtol = metric_rtol, cache = TRUE,
                            dense_maxn = maxeig, name = mat_name)
      attr(M, cache_attr) <- fac
    }
    fac$matrix <- M
    fac
  }

  # Top-k eigenpairs of a symmetric PSD target matrix. Full-rank requests (and
  # tiny problems, where dense is exact and cheap) use base::eigen.
  top_eigs <- function(target_mat, k_req, label) {
    dim_t <- nrow(target_mat)
    if (k_req >= dim_t - 1L || dim_t <= 500L) {
      es <- eigen(target_mat, symmetric = TRUE)
      keep <- seq_len(min(k_req, dim_t))
      list(values = es$values[keep], vectors = es$vectors[, keep, drop = FALSE])
    } else {
      tryCatch(.top_eigs_sym(target_mat, k_req, "LA", tol = 1e-10),
               error = function(e) stop("Eigen decomp failed in gmdLA (", label, "): ", e$message))
    }
  }

  select_valid <- function(eig_res, label) {
    valid_idx <- which(eig_res$values > 0 &
                         eig_res$values > rank_rtol^2 * max(eig_res$values, 0))
    if (length(valid_idx) == 0) stop("No positive eigenvalues found in gmdLA (", label, ").")
    list(values = eig_res$values[valid_idx],
         vectors = eig_res$vectors[, valid_idx, drop = FALSE])
  }

  # --- Main Logic --- #
  if (!use_dual) { # Primal: n_orig >= p_orig
      if (verbose) message(" gmdLA: Using primal approach (n >= p)")
      Rf <- metric_factor_cached(R, paste0(cache_attr_name, "_R"), "R")
      R <- Rf$matrix

      if (verbose) message("  Calculating X'QX...")
      XQX <- Matrix::crossprod(X, Q) %*% X # p x p matrix

      if (verbose) message("  Calculating F' X'QX F...")
      target_mat <- as.matrix(Rf$apply_t(XQX %*% Rf$mat))
      target_mat <- 0.5 * (target_mat + t(target_mat))

      if (verbose) message("  Performing eigen decomposition on target matrix (dim: ", nrow(target_mat), ")...")
      k_request <- min(k, p_orig, nrow(target_mat))
      if (k_request < 1) stop("k_request must be >= 1 in gmdLA (primal)")
      eig_res <- select_valid(top_eigs(target_mat, k_request, "primal"), "primal")
      eig_vals <- eig_res$values
      eig_vecs <- eig_res$vectors
      k_found <- length(eig_vals)
      if (verbose) message(paste("  (Found ", k_found, " eigenvalues > tol)"))

      dgmd <- sqrt(eig_vals)

      if (verbose) message("  Calculating ov (V = F^{-T} * eigenvectors)...")
      vgmd <- as.matrix(Rf$solve_t(eig_vecs)) # ov (p x k_found)

      if (verbose) message("  Calculating ou (U)...")
      # ugmd_i = X R vgmd_i / ||vgmd_i||_{RnR}, with RnR = R (X'QX) R.
      # Work with W = R vgmd (p x k) so no p x p product beyond XQX is formed:
      # the norm is (R v_i)' XQX (R v_i) and X R v_i = X W_i.
      W <- R %*% vgmd                                   # p x k_found
      norms_sq <- as.numeric(Matrix::colSums(W * (XQX %*% W)))
      ugmd <- matrix(0.0, n_orig, k_found)
      ok <- is.finite(norms_sq) & (norms_sq > rank_rtol^2 * eig_vals[1])
      if (any(!ok)) {
          warning("Near-zero norm encountered during ugmd normalization for component(s) ",
                  paste(which(!ok), collapse = ", "))
      }
      if (any(ok)) {
          XW <- as.matrix(X %*% W[, ok, drop = FALSE])  # n x sum(ok), one BLAS call
          ugmd[, ok] <- sweep(XW, 2, sqrt(norms_sq[ok]), `/`)
      }
      total_variance <- sum(XQX * R)                    # tr(X'QX R), reuses XQX
  } else { # Dual: n_orig < p_orig
      if (verbose) message(" gmdLA: Using dual approach (n < p)")
      Qf <- metric_factor_cached(Q, paste0(cache_attr_name, "_Q"), "Q")
      Q <- Qf$matrix

      if (verbose) message("  Calculating X R X'...")
      XRXt <- X %*% R %*% Matrix::t(X) # n x n matrix

      if (verbose) message("  Calculating F' X R X' F...")
      target_mat <- as.matrix(Qf$apply_t(XRXt %*% Qf$mat))
      target_mat <- 0.5 * (target_mat + t(target_mat))

      if (verbose) message("  Performing eigen decomposition on target matrix (dim: ", nrow(target_mat), ")...")
      k_request <- min(k, n_orig, nrow(target_mat))
      if (k_request < 1) stop("k_request must be >= 1 in gmdLA (dual)")
      eig_res <- select_valid(top_eigs(target_mat, k_request, "dual"), "dual")
      eig_vals <- eig_res$values
      eig_vecs <- eig_res$vectors
      k_found <- length(eig_vals)
      if (verbose) message(paste("  (Found ", k_found, " eigenvalues > tol)"))

      dgmd <- sqrt(eig_vals)

      if (verbose) message("  Calculating ou (U = F^{-T} * eigenvectors)...")
      ugmd <- as.matrix(Qf$solve_t(eig_vecs)) # ou (n x k_found)

      if (verbose) message("  Calculating ov (V = X' Q U D^{-1})...")
      # The GMD factors satisfy vgmd = X' Q ugmd D^{-1} (row metric included).
      Xt_ugmd <- Matrix::crossprod(X, Q %*% ugmd) # p x k_found
      vgmd_unnorm <- sweep(as.matrix(Xt_ugmd), 2, dgmd, `/`) # p x k_found

      # Renormalize in the R metric (numerical no-op in exact arithmetic)
      # vgmd_unnorm already has unit R-norm in exact arithmetic (X'Q u_i = d_i v_i),
      # so this guard is dimensionless: only a genuinely degenerate column fails it.
      norms_sq <- as.numeric(Matrix::colSums(vgmd_unnorm * as.matrix(R %*% vgmd_unnorm)))
      vgmd <- matrix(0.0, p_orig, k_found)
      ok <- is.finite(norms_sq) & (norms_sq > rank_rtol^2)
      if (any(!ok)) {
          warning("Near-zero norm encountered during vgmd normalization for component(s) ",
                  paste(which(!ok), collapse = ", "))
      }
      if (any(ok)) {
          vgmd[, ok] <- sweep(vgmd_unnorm[, ok, drop = FALSE], 2,
                              sqrt(norms_sq[ok]), `/`)
      }
      total_variance <- sum(Q * XRXt)             # tr(X'QX R) = tr(Q X R X'), reuses XRXt
  }

  # Calculate explained variance
  if (verbose) message(" Calculating explained variance...")
  total_variance <- as.numeric(total_variance)
  if (!is.finite(total_variance)) total_variance <- NA_real_

  if (is.na(total_variance) || total_variance <= 0) {
      propv <- rep(0, k_found)
      warning("Total generalized variance is near zero or could not be computed.")
  } else {
      propv <- (dgmd^2) / total_variance
  }
  cumv <- cumsum(propv)

  if (k_found < k) {
      warning("gmdLA: Found only ", k_found, " positive eigenvalues > tol, less than requested k=", k)
  }

  list(
    u = ugmd,       # ou (dimension depends on dual/primal)
    v = vgmd,       # ov (dimension depends on dual/primal)
    d = dgmd,       # singular values (k_found)
    k = k_found,    # Number of components found
    cumv = cumv,    # Cumulative variance explained
    propv = propv   # Proportion of variance explained
  )
}
#' @keywords internal
gmd_deflation_cpp_dispatch <- function(X, Q, R, k, thr = 1e-7, maxit = 500L,
                                       verbose = FALSE, rank_rtol = 1e-6) {
  if (methods::is(X, "sparseMatrix")) {
    if (exists("gmd_deflation_cpp_sp", mode = "function")) {
      res <- gmd_deflation_cpp_sp(as_dgc(X), as_dgc(Q), as_dgc(R), k,
                                  thr = thr, maxit = maxit, verbose = verbose, rank_rtol = rank_rtol)
      return(.emit_deflation_warnings(res))
    }
    warning("Sparse C++ deflation backend is unavailable; falling back to R deflation to avoid densifying X.")
    return(gmd_deflationR(X, Q, R, k, thr = thr, maxit = maxit, verbose = verbose, rank_rtol = rank_rtol))
  }

  res <- gmd_deflation_cpp(X, Q, R, k, thr = thr, maxit = maxit, verbose = verbose, rank_rtol = rank_rtol)
  .emit_deflation_warnings(res)
}


.emit_deflation_warnings <- function(res) {
  if (!is.null(res$warnings) && length(res$warnings)) {
    for (w in res$warnings) warning(w, call. = FALSE)
  }
  res$warnings <- NULL
  res
}


#' @noRd
#' @importFrom Matrix diag crossprod t
#' @importFrom stats rnorm
gmd_deflationR <- function(X, Q, R, k, thr = 1e-6, maxit = 500L, verbose = FALSE,
                            rank_rtol = 1e-6) {

  n <- nrow(X)
  p <- ncol(X)
  if (!is.numeric(maxit) || length(maxit) != 1 || maxit < 1 || maxit != floor(maxit)) {
    stop("maxit must be a single positive integer.")
  }
  max_iter_defl <- as.integer(maxit)

  ugmd <- matrix(0.0, n, k)
  vgmd <- matrix(0.0, p, k)
  dgmd <- numeric(k)
  propv <- numeric(k)

  # Calculate total variance once: trace(X' Q X R)
  qrnorm <- tryCatch(sum((Q %*% X) * (X %*% R)), # tr(X'QXR) without forming p x p
                    error = function(e) {
                        warning("Could not compute total variance trace: ", e$message)
                        NA_real_
                    })

  scale_ref <- 1
  if (is.na(qrnorm) || qrnorm <= 0) {
      warning("Total generalized variance is not positive or could not be computed.")
      qrnorm <- 1 # Avoid division by zero, propv will be inaccurate
  } else {
      scale_ref <- sqrt(qrnorm)
  }
  # Norm/singular-value cutoffs are relative to the scale of X in the (Q, R)
  # metric, so results are invariant to rescaling X.
  norm_floor_sq <- (.Machine$double.eps * scale_ref)^2

  k_found <- 0

  residual_mv <- function(w, n_prev) {
    y <- X %*% w
    if (n_prev > 0) {
      idx <- seq_len(n_prev)
      coeff <- dgmd[idx] * as.numeric(Matrix::crossprod(vgmd[, idx, drop = FALSE], w))
      y <- y - ugmd[, idx, drop = FALSE] %*% coeff
    }
    as.matrix(y)
  }

  residual_t_mv <- function(z, n_prev) {
    y <- Matrix::crossprod(X, z)
    if (n_prev > 0) {
      idx <- seq_len(n_prev)
      coeff <- dgmd[idx] * as.numeric(Matrix::crossprod(ugmd[, idx, drop = FALSE], z))
      y <- y - vgmd[, idx, drop = FALSE] %*% coeff
    }
    as.matrix(y)
  }

  for (i in 1:k) {
    if (verbose) message(paste(" Deflation component", i))
    # Initialize u, v for power iteration
    u <- matrix(stats::rnorm(n), ncol = 1)
    u <- u / sqrt(sum(u^2))
    v <- matrix(stats::rnorm(p), ncol = 1)
    v <- v / sqrt(sum(v^2))

    iter <- 0
    converged <- FALSE

    while (iter < max_iter_defl) {
      iter <- iter + 1
      oldu <- u
      oldv <- v

      # Update u: u_hat = X_residual R v; normalize u = u_hat / sqrt(u_hat' Q u_hat)
      uhat <- residual_mv(R %*% v, k_found)
      u_norm_sq <- as.numeric(Matrix::crossprod(uhat, Q) %*% uhat)
      if (!is.finite(u_norm_sq) || u_norm_sq <= norm_floor_sq) { # Check for near zero norm
          if (verbose) message("  u norm near zero, stopping power iteration for component ", i)
          break
      }
      u <- uhat / sqrt(as.numeric(u_norm_sq))

      # Update v: v_hat = X_residual' Q u; normalize v = v_hat / sqrt(v_hat' R v_hat)
      vhat <- residual_t_mv(Q %*% u, k_found)
      v_norm_sq <- as.numeric(Matrix::crossprod(vhat, R) %*% vhat)
       if (!is.finite(v_norm_sq) || v_norm_sq <= norm_floor_sq) { # Check for near zero norm
          if (verbose) message("  v norm near zero, stopping power iteration for component ", i)
          break
      }
      v <- vhat / sqrt(as.numeric(v_norm_sq))

      # Check convergence (squared norm difference)
      err <- sum((oldu - u)^2) + sum((oldv - v)^2)
      if (err < thr) {
          converged <- TRUE
          if (verbose) message(paste("  Power iteration converged in", iter, "steps."))
          break
      }
    } # End inner while loop

    if (!converged) {
        warning("Power iteration did not converge for component ", i, " within ", max_iter_defl, " iterations.")
        # If not converged, should we stop? Or continue with the current u,v?
        # Let's stop deflation here if power method fails
        warning("Stopping deflation due to non-convergence of power iteration.")
        break # Exit outer for loop
    }

    # Calculate singular value d = u' Q X_residual R v
    d_i <- Matrix::crossprod(u, Q) %*% residual_mv(R %*% v, k_found)
    current_d <- as.numeric(d_i)

    # Check for degenerate component: relative to the largest singular value
    # extracted so far; the first component establishes that scale.
    d_ref <- if (k_found > 0) abs(dgmd[1]) else abs(current_d)
    if (!is.finite(current_d) || current_d <= 0 || current_d <= rank_rtol * d_ref) {
        warning("Component ", i, " is degenerate (singular value near zero: ", signif(current_d, 3), "). Stopping deflation.")
        break # Exit outer for loop
    }

    # Store results for this valid component
    k_found <- k_found + 1
    dgmd[k_found] <- current_d
    ugmd[, k_found] <- u[, 1]
    vgmd[, k_found] <- v[, 1]

    # Calculate proportion of variance for this component
    propv[k_found] <- dgmd[k_found]^2 / qrnorm

  } # End outer for loop (components)

  if (k_found < k) {
      warning("Deflation stopped early. Found ", k_found, " components instead of requested ", k, ".")
      # Trim result arrays
      dgmd <- dgmd[seq_len(k_found)]
      ugmd <- ugmd[, seq_len(k_found), drop = FALSE]
      vgmd <- vgmd[, seq_len(k_found), drop = FALSE]
      propv <- propv[seq_len(k_found)]
  }

  if (k_found == 0) {
      return(list(d = numeric(0),
                  v = matrix(0, p, 0),
                  u = matrix(0, n, 0),
                  k = 0,
                  cumv = numeric(0),
                  propv = numeric(0)))
  }

  cumv <- cumsum(propv) # Calculate cumulative sum on valid components

  list(d = as.vector(dgmd), v = vgmd, u = ugmd, k = k_found, cumv = cumv, propv = propv)
}


#' Truncate a genpca fit to fewer components
#'
#' Returns a new `genpca` object retaining only the first `ncomp` components.
#' All component-indexed slots (`v`, `s`, `sdev`, `ov`, `ou`, `u`, `propv`,
#' `cumv`) are sliced consistently; the preprocessing object and constraint
#' matrices are carried over unchanged.
#'
#' @param x A `genpca` object.
#' @param ncomp Number of components to retain (a positive integer no larger
#'   than `ncomp(x)`).
#' @return A `genpca` object with `ncomp` components.
#' @seealso [genpca()], [reconstruct.genpca()]
#' @examples
#' X <- matrix(rnorm(60), 15, 4)
#' fit <- genpca(X, ncomp = 4)
#' fit2 <- truncate(fit, 2)
#' multivarious::ncomp(fit2)
#' @importFrom multivarious ncomp scores sdev bi_projector
#' @export
truncate.genpca <- function(x, ncomp) {
  # Check requested ncomp
  current_ncomp <- multivarious::ncomp(x)
  if (missing(ncomp)) stop("Argument 'ncomp' must be provided.")
  if (!is.numeric(ncomp) || length(ncomp) != 1 || ncomp < 1 || ncomp > current_ncomp || ncomp != floor(ncomp)) {
      stop(paste0("Requested ncomp (", ncomp, ") must be a positive integer <= ", current_ncomp, "."))
  }

  if (ncomp == current_ncomp) return(x) # Nothing to do

  # Use the bi_projector constructor to create the truncated object
  # Select the first 'ncomp' components from relevant slots
  ret <- multivarious::bi_projector(
    v = x$v[, 1:ncomp, drop = FALSE], # A ov
    s = multivarious::scores(x)[, 1:ncomp, drop = FALSE],   # ou D
    sdev = multivarious::sdev(x)[1:ncomp],                # d
    preproc = x$preproc,                    # Preprocessing object
    ov = x$ov[, 1:ncomp, drop = FALSE],       # Orthonormal V
    ou = x$ou[, 1:ncomp, drop = FALSE],       # Orthonormal U
    u = x$u[, 1:ncomp, drop = FALSE],         # M ou
    classes = "genpca",                     # bi_projector() appends "bi_projector", "projector"
    A = x$A,                                # Constraint matrix A
    M = x$M,                                # Constraint matrix M
    propv = if (!is.null(x$propv)) x$propv[1:ncomp] else NULL, # Proportion variance
    cumv = if (!is.null(x$cumv)) x$cumv[1:ncomp] else NULL    # Cumulative variance
  )
  return(ret) # Explicitly return the new object
}


#' Reconstruct data from a genpca fit
#'
#' Reconstructs (an approximation of) the original data from a `genpca` fit as
#' `ou[, comp] %*% diag(d[comp]) %*% t(ov[, comp])`, followed by the inverse of
#' the preprocessing transform. With all components and full rank this
#' recovers the original data.
#'
#' @param x A `genpca` object.
#' @param comp Integer vector of components to use (default: all).
#' @param rowind Optional integer vector of rows to reconstruct (default: all).
#' @param colind Optional integer vector of columns to reconstruct (default:
#'   all). The inverse preprocessing transform is applied to the selected
#'   columns.
#' @param ... Ignored.
#' @return A numeric matrix of dimension `length(rowind) x length(colind)`.
#' @seealso [genpca()], [truncate.genpca()]
#' @examples
#' X <- matrix(rnorm(60), 15, 4)
#' fit <- genpca(X, ncomp = 4, preproc = multivarious::center())
#' max(abs(reconstruct(fit) - X)) # ~ 0 at full rank
#' @importFrom multivarious ncomp sdev scores inverse_transform
#' @importFrom assertthat assert_that
#' @importFrom Matrix Diagonal t
#' @export
reconstruct.genpca <- function(x,
                               comp = 1:multivarious::ncomp(x),
                               rowind = NULL, # Default to all rows
                               colind = NULL, # Default to all cols
                               ...) {

  max_comp <- multivarious::ncomp(x)
  if (max_comp == 0) return(matrix(0,
                                   nrow = length(rowind %||% 1:nrow(x$M)),
                                   ncol = length(colind %||% 1:nrow(x$A)))) # Return empty if no components
  if (length(comp) == 0 || min(comp) < 1 || max(comp) > max_comp) {
      stop("Selected components 'comp' are out of bounds [1, ", max_comp, "].")
  }

  # Reconstruction uses U D V' where U,V are orthonormal in M,A metrics (ou, ov)
  # X_hat_preproc = ou[, comp] %*% D[comp, comp] %*% t(ov[, comp])

  dvals <- multivarious::sdev(x)[comp]
  # Use Matrix::Diagonal for efficiency
  D_comp <- Matrix::Diagonal(n = length(dvals), x = dvals)

  # Determine effective row/col indices for ou/ov matrices
  eff_rowind <- rowind %||% 1:nrow(x$ou)
  eff_colind <- colind %||% 1:nrow(x$ov)

  # Helper function for safe indexing
  safe_index <- function(mat, rows, cols) {
     if (is.null(rows) && is.null(cols)) return(mat)
     if (is.null(rows)) rows <- 1:nrow(mat)
     if (is.null(cols)) cols <- 1:ncol(mat)
     mat[rows, cols, drop = FALSE]
  }

  # Select the specified components and indices for ou and ov
  OU_comp <- safe_index(x$ou, eff_rowind, comp)
  OV_comp <- safe_index(x$ov, eff_colind, comp) # ov rows correspond to X columns

  # Perform the core reconstruction: OU %*% D %*% t(OV)
  # Ensure matrix multiplication handles sparse matrices correctly
  reconstructed_data_preproc <- OU_comp %*% D_comp %*% Matrix::t(OV_comp)

  # Convert to regular matrix for inverse_transform (it requires a matrix, not Matrix object)
  reconstructed_data_preproc <- as.matrix(reconstructed_data_preproc)

  # Apply inverse pre-processing transform. colind must be forwarded so that
  # e.g. centering adds back the means of the *selected* columns, not the
  # first length(colind) ones.
  final_reconstruction <- multivarious::inverse_transform(
    x$preproc, reconstructed_data_preproc, colind = colind)

  return(final_reconstruction)
}

# Helper for default NULL indexing
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}


# S3 method for ncomp (registered via S3method; not exported by name)
#' @noRd
#' @export
ncomp.genpca <- function(x) {
  # Number of components is determined by the length of singular values
  # or columns in ou/ov/s/v
  length(multivarious::sdev(x))
}
