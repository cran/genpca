# Import necessary functions from packages
#' @importFrom Matrix Matrix Diagonal crossprod tcrossprod t solve bandSparse sparseMatrix Cholesky rowSums diag<-
#' @importFrom FNN get.knn
#' @importFrom stats median
NULL

# Function to construct the second differences matrix
#' @noRd
second_diff_matrix <- function(n) {
  # n: Length of the time series
  if (n < 3) stop("n must be at least 3 to compute second differences.")

  # Number of second differences is (n - 2)
  m <- n - 2

  # Create the diagonals for the second difference matrix
  # Diagonals for D2: (1, -2, 1) pattern for n-2 rows, n cols
  # k=0: 1s at (1,1), (2,2), ..., (m,m)
  # k=1: -2s at (1,2), (2,3), ..., (m,m+1)
  # k=2: 1s at (1,3), (2,4), ..., (m,n)
  diagonals <- list(
    rep(1, m),   # for k=0 offset
    rep(-2, m),  # for k=1 offset
    rep(1, m)    # for k=2 offset
  )

  # Build the (n - 2) x n second difference matrix using bandSparse
  # Correct offsets k = c(0, 1, 2) relative to main diagonal
  D2 <- Matrix::bandSparse(m, n, k = 0:2, diagonals = diagonals, symmetric = FALSE)

  return(D2)
}




#' Sparse and Functional Principal Components Analysis (SFPCA) with Spatial Coordinates
#'
#' Performs Sparse and Functional PCA on a data matrix, allowing for both sparsity and smoothness
#' in the estimated principal components. Penalty parameters left `NULL` are selected
#' automatically (see Details). The spatial smoothness penalty is constructed based on
#' provided spatial coordinates.
#'
#' Each rank-1 problem is solved by alternating solves of the penalized
#' quadratic subproblems (via C++ coordinate descent) followed by rescaling
#' onto the smoothness-metric ball, in the constraint form of Allen & Weylandt
#' (2019). For the convex `"l1"` penalty with subproblems solved to tolerance
#' (the internal `exact_inner = TRUE` path, used by the monotonicity test) the
#' objective is monotonically non-decreasing; the default inexact path
#' tightens the inner tolerance to a floor before it may declare convergence,
#' reproducing the same terminal iterates but without an every-iteration
#' monotonicity guarantee (it may also stop at `max_iter`).
#'
#' When `lambda_u` or `lambda_v` is `NULL` it is selected per component by a
#' BIC-style criterion along a regularization path. For the convex `"l1"`
#' penalty `lambda_max = max(abs(b))` is, in closed form, the smallest value
#' whose subproblem solution is exactly zero (at `x = 0` the `S x` term
#' vanishes, so the KKT condition `|b_j| <= lambda` does not depend on `S`);
#' where `b` is the matrix-vector product with the other factor fixed at the
#' SVD initializer. For the non-convex `"scad"` penalty the same value anchors
#' the path but is not a global-optimality threshold. `nlambda` values are
#' laid log-spaced down to `lambda_min_ratio * lambda_max`, coordinate descent
#' is warm-started along the path, and the value minimizing
#' `log(RSS / (n p)) + df * log(n p) / (n p)` is chosen, with `df` the support
#' size of the solution and `RSS` the one-sided rank-1 residual sum of squares
#' with the opposite factor held fixed (a selection heuristic, not the BIC of
#' the fully alternated rank-1 model). The all-zero solution (at `lambda_max`)
#' is a legitimate candidate: if no rank-1 structure justifies its degrees of
#' freedom, the component is returned as exactly zero with `d = 0`.
#'
#' When `alpha_u` or `alpha_v` is `NULL` it defaults to
#' `1 / lambda_max(Omega)`, so the roughest direction of the smoothness
#' penalty is weighted exactly as strongly as the identity term. This makes
#' the default invariant to the scaling of `Omega` and bounds the condition
#' number of every subproblem system `I + alpha * Omega` by 2.
#'
#' @param X A numeric data matrix of dimensions n (observations/time points) by p (variables/space).
#' @param K The number of principal components to estimate.
#' @param spat_cds A matrix of spatial coordinates for each column of X (variables). Each row
#'   corresponds to a spatial dimension (e.g., x, y, z), and each column corresponds to a variable.
#'   Note the orientation: this is `dimensions x variables`, so
#'   `ncol(spat_cds)` must equal `ncol(X)` -- the transpose of the layout a
#'   coordinate data frame usually has. For a one-dimensional axis (a spectrum,
#'   a transect) pass `matrix(coords, nrow = 1)`.
#' @param lambda_u Sparsity penalty parameter for u. If NULL, selected per component by BIC
#'   along a regularization path (see Details).
#' @param lambda_v Sparsity penalty parameter for v. If NULL, selected per component by BIC
#'   along a regularization path (see Details).
#' @param alpha_u Smoothness penalty parameter for u. If NULL, defaults to
#'   `1 / lambda_max(Omega_u)` (see Details).
#' @param alpha_v Smoothness penalty parameter for v. If NULL, defaults to
#'   `1 / lambda_max(Omega_v)` (see Details).
#' @param Omega_u A positive semi-definite matrix for smoothness penalty on u. If NULL, defaults to
#'   second differences penalty (sparse matrix). Unlike `Omega_u`, there is no
#'   corresponding `Omega_v` argument: the column-side smoothness penalty is
#'   always built internally from `spat_cds` (via `knn`); supplying a custom
#'   `Omega_v` is not currently supported.
#' @param penalty_u The penalty function for u. Either "l1" (lasso, the
#'   default) or "scad".
#' @param penalty_v The penalty function for v. Either "l1" (lasso, the
#'   default) or "scad".
#' @param nlambda Number of values on the regularization path used for BIC selection
#'   of `lambda_u`/`lambda_v` when they are NULL. Default `10`.
#' @param lambda_min_ratio Smallest path value as a fraction of the closed-form
#'   `lambda_max`, on a log-spaced grid. Default `1e-2`.
#' @param knn Number of nearest neighbours for constructing `Omega_v`. Default
#'   `min(6, ncol(X) - 1)`.
#' @param max_iter Maximum number of iterations for the alternating optimization.
#'   Default `100`.
#' @param tol Tolerance for convergence of the rank-1 objective. Default `1e-6`.
#' @param verbose Logical; if TRUE, prints progress messages.
#' @param uthresh Deprecated and ignored; `lambda_u` is now selected by BIC.
#' @param vthresh Deprecated and ignored; `lambda_v` is now selected by BIC.
#' @return An object of class `c("sfpca", "bi_projector")` from the
#'   \pkg{multivarious} framework. Use `multivarious::scores()` for the sample
#'   scores (\eqn{U D}), `multivarious::components()` for the sparse loadings
#'   \eqn{V}, `multivarious::sdev()` for \eqn{d_k}, and
#'   `multivarious::reconstruct()` for the rank-`K` approximation. `ov` (like
#'   `components()`) holds the sparse right factors \eqn{V}; `ou` holds the
#'   left factors \eqn{U}. The selected penalty parameters are stored as
#'   `lambda_u`, `lambda_v`, `alpha_u`, and `alpha_v`. For backward
#'   compatibility the pre-0.1 list fields `$d` (singular values) and `$u`
#'   (left factors) remain readable but emit a deprecation warning; use
#'   `sdev()` and `scores()`/`$ou` instead.
#'
#'   \strong{Important:} unlike `genpca()`, the columns of `U` (`ou`) and `V`
#'   (`ov`) are Euclidean unit-norm but are **not** mutually orthogonal
#'   across components -- `sfpca()` extracts each rank-1 term from a
#'   constraint-form subproblem rather than a joint SVD, so `U'U != I` and
#'   `V'V != I` in general. Consequently `multivarious::sdev()` here is
#'   *not* the singular values of `X`; it is the per-component captured
#'   covariance \eqn{d_k = u_k' X_k v_k}, where \eqn{X_k} is the matrix after
#'   the preceding components have been deflated out (so the identity holds
#'   against `X` itself only for \eqn{k = 1}). This non-orthogonality is also why
#'   `reconstruct()` for `"sfpca"` objects uses the stored `U`, `d`, `V`
#'   factors directly (`U D V'`) rather than SVD-based identities such as the
#'   Moore-Penrose pseudoinverse of the loadings, which would not reproduce
#'   the fitted model for non-orthogonal `V` (see `reconstruct.sfpca()`).
#' @references Allen, G. I., & Weylandt, M. (2019). Sparse and functional
#'   principal components analysis. In \emph{2019 IEEE Data Science Workshop
#'   (DSW)} (pp. 11-16). \doi{10.1109/DSW.2019.8755778}. Also available as
#'   \href{https://arxiv.org/abs/1309.2895}{arXiv:1309.2895}, first posted in
#'   2013 and revised through 2019; the preprint and the DSW paper are the
#'   same work, which is why both years appear in the literature.
#' @seealso [genpca()] for the shared \pkg{multivarious} verbs;
#'   \code{multivarious::bi_projector}.
#' @examples
#' library(Matrix)
#' set.seed(123)
#' # Smooth temporal factor, sparse spatial factor
#' n <- 100  # Number of time points
#' p <- 50   # Number of spatial locations
#' u <- sin(seq(0, 2 * pi, length.out = n))
#' v <- c(rnorm(10), rep(0, p - 10))
#' X <- 8 * tcrossprod(u / sqrt(sum(u^2)), v / sqrt(sum(v^2))) +
#'   matrix(rnorm(n * p, sd = 0.2), n, p)
#' spat_cds <- matrix(runif(p * 3), nrow = 3, ncol = p)  # 3D coordinates
#' result <- sfpca(X, K = 1, spat_cds = spat_cds)
#' multivarious::sdev(result)                  # captured covariance (BIC-tuned)
#' sum(multivarious::components(result) != 0)  # sparse spatial loading
#' @export
sfpca <- function(X, K, spat_cds,
                  lambda_u = NULL, lambda_v = NULL,
                  alpha_u = NULL, alpha_v = NULL,
                  Omega_u = NULL,
                  penalty_u = "l1", penalty_v = "l1",
                  nlambda = 10, lambda_min_ratio = 1e-2,
                  knn = min(6, ncol(X) - 1),  # Number of nearest neighbors for spatial penalty
                  max_iter = 100, tol = 1e-6, verbose = FALSE,
                  uthresh = NULL, vthresh = NULL) {

  if (!is.null(uthresh) || !is.null(vthresh)) {
    warning("`uthresh` and `vthresh` are deprecated and ignored: ",
            "`lambda_u`/`lambda_v` are now selected by BIC along a ",
            "regularization path. Pass them explicitly to fix the penalties.")
  }
  if (length(nlambda) != 1 || !is.finite(nlambda) || nlambda < 2) {
    stop("`nlambda` must be a single number >= 2.")
  }
  nlambda <- as.integer(nlambda)
  if (length(lambda_min_ratio) != 1 || !is.finite(lambda_min_ratio) ||
      lambda_min_ratio <= 0 || lambda_min_ratio >= 1) {
    stop("`lambda_min_ratio` must be in (0, 1).")
  }

  n <- nrow(X)
  p <- ncol(X)
  d_list <- numeric(K)
  u_list <- matrix(0, n, K)
  v_list <- matrix(0, p, K)
  lambda_u_list <- numeric(K)
  lambda_v_list <- numeric(K)
  alpha_u_list <- numeric(K)
  alpha_v_list <- numeric(K)

  # Convert X to Matrix if not already
 if (!inherits(X, "Matrix")) {
    X <- Matrix::Matrix(X, sparse = FALSE)
  }

  # Default Omega_u if not provided (second differences)
  if (is.null(Omega_u)) {
    D_n <- second_diff_matrix(n)  # Already returns sparse Matrix
    Omega_u <- Matrix::crossprod(D_n)
  }

  # Construct Omega_v based on spatial coordinates
  # We will use a weighted graph Laplacian where weights are based on distances
  if (is.null(spat_cds)) {
    stop("spat_cds must be provided for constructing the spatial penalty matrix Omega_v.")
  }

  # `spat_cds` is dimensions x variables, which is the transpose of the layout
  # most callers reach for. Catch the mismatch here: downstream it silently
  # redefines p and fails much later with an unrelated message.
  if (is.null(dim(spat_cds)) || length(dim(spat_cds)) != 2L) {
    stop("`spat_cds` must be a matrix with one column per variable ",
         "(spatial dimensions in rows, variables in columns). For a ",
         "one-dimensional axis use matrix(coords, nrow = 1).")
  }
  if (ncol(spat_cds) != p) {
    msg <- paste0("`spat_cds` must have one column per variable: ncol(X) is ",
                  p, " but ncol(spat_cds) is ", ncol(spat_cds), ".")
    if (nrow(spat_cds) == p) {
      msg <- paste0(msg, " It looks transposed -- `spat_cds` is ",
                    "dimensions x variables, so pass t(spat_cds).")
    }
    stop(msg)
  }
  if (anyNA(spat_cds)) {
    stop("`spat_cds` must not contain NA: nearest-neighbour distances are undefined.")
  }

  if (verbose) cat("Constructing spatial penalty matrix Omega_v based on spat_cds...\n")

  Omega_v <- construct_spatial_penalty(spat_cds, k = knn)

  # Smoothness weights: scale-free default alpha = 1 / lambda_max(Omega), so
  # the roughest penalty direction is weighted exactly as strongly as the
  # identity term and cond(I + alpha * Omega) <= 2 however Omega was scaled.
  alpha_u_use <- if (is.null(alpha_u)) default_alpha(Omega_u) else alpha_u
  alpha_v_use <- if (is.null(alpha_v)) default_alpha(Omega_v) else alpha_v
  # S = I + alpha * Omega must stay SPD for the coordinate-descent subproblems
  # to be convex; a negative or non-finite alpha breaks that precondition (the
  # C++ solver only checks a positive diagonal, which is not sufficient).
  for (nm in c("alpha_u", "alpha_v")) {
    a <- get(paste0(nm, "_use"))
    if (length(a) != 1 || !is.finite(a) || a < 0) {
      stop("`", nm, "` must be a single finite non-negative number.")
    }
  }
  if (verbose && (is.null(alpha_u) || is.null(alpha_v))) {
    cat("Default smoothness weights: alpha_u =", alpha_u_use,
        "alpha_v =", alpha_v_use, "\n")
  }
  S_u <- as_dgc(Matrix::Diagonal(n) + alpha_u_use * Omega_u)
  S_v <- as_dgc(Matrix::Diagonal(p) + alpha_v_use * Omega_v)

  # BIC selection of the sparsity penalties needs ||X_res||_F^2 per component.
  need_lambda_path <- is.null(lambda_u) || is.null(lambda_v)
  F2_X <- if (need_lambda_path) sum(X^2) else NA_real_

  # Implicit deflation: X is never modified (sparse inputs stay sparse).
  # Extracted components deflate the operator products instead:
  # X_res v = Xv - U diag(d) (V'v).
  for (k in 1:K) {
    if (verbose) cat("Component", k, "\n")
    U_prev <- if (k > 1) u_list[, seq_len(k - 1), drop = FALSE] else NULL
    V_prev <- if (k > 1) v_list[, seq_len(k - 1), drop = FALSE] else NULL
    d_prev <- if (k > 1) d_list[seq_len(k - 1)] else numeric(0)
    ops <- sfpca_make_ops(X, U_prev, d_prev, V_prev)
    # Estimate initial u and v using SVD of the (implicitly) deflated matrix
    svd_res <- svd1_deflated(X, U_prev, d_prev, V_prev, verbose = verbose)
    u_init <- as.numeric(svd_res$u)
    v_init <- as.numeric(svd_res$v)
    d_init <- svd_res$d[1]
    # Sparsity penalties: per-component BIC along a warm-started lambda path.
    # The fixed factor is normalized exactly as the first alternation step
    # normalizes it (S-norm), so the selected lambda lives on the same scale
    # as the subproblems it will be used in.
    if (need_lambda_path) {
      F2_res <- sfpca_res_fnorm2(F2_X, X, U_prev, d_prev, V_prev)
    }
    if (is.null(lambda_u)) {
      v_fix <- v_init / s_norm(v_init, S_v)
      sel_u <- sfpca_select_lambda(ops$mv(v_fix), S_u, F2_res, n * p,
                                   sqrt(sum(v_fix^2)), penalty_u,
                                   nlambda, lambda_min_ratio)
      lambda_u_k <- sel_u$lambda
    } else {
      lambda_u_k <- lambda_u
    }
    if (is.null(lambda_v)) {
      u_fix <- u_init / s_norm(u_init, S_u)
      sel_v <- sfpca_select_lambda(ops$tmv(u_fix), S_v, F2_res, n * p,
                                   sqrt(sum(u_fix^2)), penalty_v,
                                   nlambda, lambda_min_ratio)
      lambda_v_k <- sel_v$lambda
    } else {
      lambda_v_k <- lambda_v
    }

    # Store penalty parameters
    lambda_u_list[k] <- lambda_u_k
    lambda_v_list[k] <- lambda_v_k
    alpha_u_list[k] <- alpha_u_use
    alpha_v_list[k] <- alpha_v_use

    if (verbose) {
      cat("Penalties for component", k, ":\n")
      cat("lambda_u =", lambda_u_k, "lambda_v =", lambda_v_k, "\n")
      cat("alpha_u =", alpha_u_use, "alpha_v =", alpha_v_use, "\n")
    }

    # Run rank-1 SFPCA with the selected penalties
    result <- sfpca_rank1(ops = ops,
                          lambda_u = lambda_u_k, lambda_v = lambda_v_k,
                          alpha_u = alpha_u_use, alpha_v = alpha_v_use,
                          Omega_u = Omega_u, Omega_v = Omega_v,
                          penalty_u = penalty_u, penalty_v = penalty_v,
                          max_iter = max_iter, tol = tol, verbose = verbose,
                          u_init = u_init, v_init = v_init, d_init = d_init)


    d_list[k] <- result$d
    u_list[, k] <- as.vector(result$u)
    v_list[, k] <- as.vector(result$v)
  }

  comp_names <- paste0("PC", seq_len(K))
  rn <- if (!is.null(rownames(X))) rownames(X) else paste0("Obs", seq_len(n))
  cn <- if (!is.null(colnames(X))) colnames(X) else paste0("Var", seq_len(p))
  dimnames(u_list) <- list(rn, comp_names)
  dimnames(v_list) <- list(cn, comp_names)

  # Scores of the sfpca deflation model X ~= U D V': F = U D (samples in
  # component space), exposed via multivarious `scores()`; loadings are the
  # sparse right factors V (`components()`). The columns of `u_list`/`v_list`
  # are Euclidean unit vectors but, unlike an SVD, are NOT mutually orthogonal
  # across components, so reconstruction uses t(V) (see reconstruct.sfpca),
  # not the pseudoinverse. sweep avoids the diag(d) scalar trap when K == 1.
  scores_mat <- sweep(u_list, 2, d_list, `*`)

  # sfpca does no preprocessing; build a fitted pass() pre_processor the same
  # way genpca does (avoids multivarious's deprecated prep() default).
  procres <- multivarious::fit_transform(multivarious::pass(),
                                         matrix(0, nrow = 1L, ncol = p))$preproc

  multivarious::bi_projector(
    v = v_list,             # loadings / components (right factors, unit-norm)
    s = scores_mat,         # scores = U D
    sdev = d_list,          # singular values
    preproc = procres,      # identity preprocessing
    ou = u_list,            # left factors (unit-norm columns, M = I here)
    ov = v_list,            # right factors (unit-norm columns, A = I here)
    lambda_u = lambda_u_list, lambda_v = lambda_v_list,
    alpha_u = alpha_u_list, alpha_v = alpha_v_list,
    classes = "sfpca"
  )
}

#' Reconstruct data from an sfpca fit
#'
#' Reconstructs the rank-`K` \code{\link{sfpca}} model as \eqn{U D V'},
#' using the stored (non-orthogonal) factors directly.
#'
#' @details
#' sfpca components are Euclidean unit vectors but are **not** mutually
#' orthogonal, so \eqn{V'V \ne I}. The inherited `reconstruct.bi_projector()`
#' method reconstructs through the Moore-Penrose pseudoinverse of the
#' loadings (`scores \%*\% pinv(V)`), which for non-orthogonal `V` does
#' **not** return the rank-`comp` model \eqn{U D V'} that `sfpca()` actually
#' fits and deflates with. This method instead computes
#' `scores(x)[rowind, comp] \%*\% t(components(x)[colind, comp])`, i.e.
#' \eqn{U D V'} restricted to the requested rows/columns/components.
#' `sfpca()` does no preprocessing, so no inverse transform is applied.
#'
#' @param x An `sfpca` object.
#' @param comp Integer vector of components to use (default: all).
#' @param rowind Optional integer vector of rows to reconstruct (default: all).
#' @param colind Optional integer vector of columns to reconstruct (default:
#'   all).
#' @param ... Ignored.
#' @return A numeric matrix of dimension `length(rowind) x length(colind)`,
#'   the rank-`length(comp)` reconstruction \eqn{U D V'} using sfpca's
#'   stored non-orthogonal factors.
#' @seealso [sfpca()]
#' @exportS3Method
reconstruct.sfpca <- function(x, comp = seq_len(multivarious::ncomp(x)),
                              rowind = NULL,
                              colind = NULL,
                              ...) {
  ncomp <- multivarious::ncomp(x)
  if (any(comp < 1) || any(comp > ncomp)) {
    stop("`comp` must index existing components (1:", ncomp, ").")
  }
  S_all <- multivarious::scores(x)
  V_all <- multivarious::components(x)
  rowind <- if (is.null(rowind)) seq_len(nrow(S_all)) else rowind
  colind <- if (is.null(colind)) seq_len(nrow(V_all)) else colind
  S <- S_all[rowind, comp, drop = FALSE]      # U D
  V <- V_all[colind, comp, drop = FALSE]      # loadings
  S %*% t(V)
}

# Deprecation shim for the pre-0.1 list return shape. `sfpca()` now returns a
# `multivarious::bi_projector`; the singular values live under `sdev()` and the
# left factors under `scores()`/`x$ou`. `$d` and `$u` remain readable but warn.
#' @noRd
.sfpca_deprecate_field <- function(field, alt) {
  warning("`$", field, "` on an sfpca object is deprecated; use ", alt, ".",
          call. = FALSE)
}

#' @rawNamespace S3method("$", sfpca)
`$.sfpca` <- function(x, name) {
  if (identical(name, "d")) {
    .sfpca_deprecate_field("d", "sdev(x)")
    return(.subset2(x, "sdev"))
  }
  if (identical(name, "u")) {
    .sfpca_deprecate_field("u", "x$ou (left factors) or scores(x) for U D")
    return(.subset2(x, "ou"))
  }
  .subset2(x, name)
}

#' Print an sfpca fit
#'
#' Prints a one-line summary of an \code{\link{sfpca}} object: number of
#' components, dimensions, and singular values.
#'
#' @param x An `sfpca` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @exportS3Method
print.sfpca <- function(x, ...) {
  k <- length(.subset2(x, "sdev"))
  cat("Sparse Functional PCA (sfpca)\n")
  cat("  components:", k, "\n")
  cat("  dims: ", nrow(.subset2(x, "ou")), " obs x ",
      nrow(.subset2(x, "v")), " vars\n", sep = "")
  cat("  singular values:",
      paste(signif(.subset2(x, "sdev"), 4), collapse = ", "), "\n")
  cat("  verbs: scores(), components(), sdev(), reconstruct()\n")
  invisible(x)
}

# Function to construct the spatial penalty matrix Omega_v based on spat_cds
#' @noRd
construct_spatial_penalty <- function(spat_cds, method = "distance", k = 6L) {
  p <- ncol(spat_cds)
  if (p <= k) {
     warning(paste0("Number of spatial locations (p=", p, ") is less than or equal to knn (k=", k, "). ",
                    "Reducing k to p-1 = ", p - 1, "."))
     k <- p - 1
  }
  if (k < 1) {
      stop("knn (k) must be at least 1.")
  }

  if (method == "distance") {
    # Compute pairwise distances between spatial coordinates is too slow for large p
    # Build adjacency matrix based on k-nearest neighbors
    knn_res <- FNN::get.knn(t(spat_cds), k = k)
    indices <- knn_res$nn.index
    distances_knn <- knn_res$nn.dist

    # Heat-kernel weights exp(-d^2 / (2 sigma^2)) with sigma set to the median
    # kNN distance. Bounded in (0, 1], so the Laplacian (and hence
    # S_v = I + alpha * Omega_v) stays well-conditioned; unbounded inverse
    # distance weights make the penalized subproblems numerically unsolvable
    # when spatial points nearly coincide.
    sigma <- stats::median(distances_knn)
    if (!is.finite(sigma) || sigma <= 0) sigma <- 1
    weights <- exp(-as.vector(t(distances_knn))^2 / (2 * sigma^2))

    # Create sparse matrix W with weights
    W <- Matrix::sparseMatrix(i = rep(1:p, each = k),
                              j = as.vector(t(indices)),
                              x = weights,
                              dims = c(p, p), index1 = TRUE) # Ensure index1=TRUE

    # Make W symmetric: W = (W + W^T) / 2
    # This ensures that if j is a neighbor of i, i is also considered a neighbor of j
    # with potentially averaged weight.
    W <- (W + Matrix::t(W)) / 2
    # Ensure diagonal is zero after symmetrization
    Matrix::diag(W) <- 0

  } else if (method == "neighbor") {
    # If a neighbor graph is provided, construct W based on it
    # This part can be customized based on available neighbor information
    stop("Neighbor graph method not implemented in this example.")
  } else {
    stop("Invalid method for constructing spatial penalty.")
  }

  # Construct graph Laplacian L = D - W
  # D is the diagonal matrix of degrees (row sums of W)
  D <- Matrix::Diagonal(x = Matrix::rowSums(W))
  L <- D - W

  # Omega_v is the graph Laplacian
  Omega_v <- L

  return(Omega_v)
}

# Matrix-vector operators for the implicitly deflated matrix
# X_res = X - U diag(d) V'. X itself is never modified, so sparse inputs
# stay sparse; each product costs O(nnz(X) + (n + p) k).
#' @noRd
sfpca_make_ops <- function(X, U = NULL, d = numeric(0), V = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  if (is.null(U) || length(d) == 0) {
    list(n = n, p = p,
         mv = function(v, args = NULL) as.numeric(X %*% v),
         tmv = function(u, args = NULL) as.numeric(Matrix::crossprod(X, u)))
  } else {
    list(n = n, p = p,
         mv = function(v, args = NULL) {
           as.numeric(X %*% v) - as.numeric(U %*% (d * crossprod(V, v)))
         },
         tmv = function(u, args = NULL) {
           as.numeric(Matrix::crossprod(X, u)) -
             as.numeric(V %*% (d * crossprod(U, u)))
         })
  }
}

# Leading singular triplet of the implicitly deflated matrix, with a dense
# fallback when the iterative solver fails
#' @noRd
svd1_deflated <- function(X, U = NULL, d = numeric(0), V = NULL,
                          verbose = FALSE) {
  ops <- sfpca_make_ops(X, U, d, V)
  tryCatch({
    sv <- .top_svd(ops$mv, 1, nu = 1, nv = 1, adjoint = ops$tmv, dim = c(ops$n, ops$p))
    if (!isTRUE(sv$converged)) stop("iterative svd did not converge")
    sv
  }, error = function(e) {
    if (verbose) cat("iterative svd failed, falling back to base R svd\n")
    Xd <- as.matrix(X)
    if (!is.null(U) && length(d) > 0) {
      Xd <- Xd - as.matrix(U) %*% (d * t(as.matrix(V)))
    }
    base_svd <- svd(Xd, nu = 1, nv = 1)
    list(u = base_svd$u[, 1, drop = FALSE],
         v = base_svd$v[, 1, drop = FALSE],
         d = base_svd$d[1])
  })
}

# Scale-free default smoothness weight alpha = 1 / lambda_max(Omega): the
# roughest direction of the penalty gets the same weight as the identity
# term, so cond(I + alpha * Omega) <= 2 no matter how Omega was scaled.
#' @noRd
default_alpha <- function(Omega) {
  lam <- tryCatch({
    if (nrow(Omega) <= 200) {
      max(eigen(as.matrix(Omega), symmetric = TRUE, only.values = TRUE)$values)
    } else {
      as.numeric(.top_eigs_sym(as_dgc(Omega), 1, "LM")$values[1])
    }
  }, error = function(e) NA_real_)
  if (!is.finite(lam)) {
    # Gershgorin upper bound for symmetric Omega (conservative: smaller alpha)
    lam <- max(Matrix::rowSums(abs(Omega)))
  }
  if (lam > 0) 1 / lam else 0
}

# Squared Frobenius norm of the implicitly deflated matrix,
# ||X - U diag(d) V'||_F^2, without forming the residual.
#' @noRd
sfpca_res_fnorm2 <- function(F2_X, X, U = NULL, d = numeric(0), V = NULL) {
  if (is.null(U) || length(d) == 0) return(F2_X)
  XV <- as.matrix(X %*% V)
  cross <- sum(d * colSums(as.matrix(U) * XV))
  Gu <- crossprod(as.matrix(U))
  Gv <- crossprod(as.matrix(V))
  # The three terms can be large and nearly cancel at high signal scale; the
  # true residual norm is non-negative, so clamp away catastrophic-cancellation
  # noise before it reaches log(RSS) downstream.
  max(F2_X - 2 * cross + sum(outer(d, d) * Gu * Gv), 0)
}

# Sparsity penalty selection for one subproblem side. lambda_max = ||b||_inf
# is the smallest lambda whose solution is exactly zero (the KKT condition at
# x = 0 does not involve S), a log-spaced path is laid down from it,
# coordinate descent is warm-started along the path, and lambda is picked by
# BIC on the rank-1 fit, log(RSS / np) + df * log(np) / np with df = support
# size (Allen & Weylandt, 2019). `fixed_norm` is the Euclidean norm of the
# fixed factor that produced b = X_res v (or X_res' u), so the fitted
# singular value is computed with unit-norm factors even though the fixed
# factor is S-normalized to match the alternation loop's scaling.
#' @noRd
sfpca_select_lambda <- function(b, S, F2, np, fixed_norm, penalty,
                                nlambda = 10, lambda_min_ratio = 1e-2) {
  lam_max <- max(abs(b))
  if (!is.finite(lam_max) || lam_max <= 0) {
    return(list(lambda = 0, lambdas = numeric(0), bic = numeric(0),
                index = NA_integer_))
  }
  lambdas <- exp(seq(log(lam_max), log(lam_max * lambda_min_ratio),
                     length.out = nlambda))
  bic <- numeric(nlambda)
  x_warm <- numeric(length(b))
  for (i in seq_len(nlambda)) {
    # Selection only ranks candidates: looser tolerance than the final solves
    sol <- sfpca_cd_solve(S, b, x_warm, lambdas[i], penalty,
                          max_sweeps = 250L, tol = 1e-5)
    x_warm <- sol$x
    df <- sum(sol$x != 0)
    if (df == 0) {
      rss <- F2
    } else {
      xn <- sol$x / sqrt(sum(sol$x^2))
      d_fit <- sum(xn * b) / fixed_norm
      rss <- max(F2 - d_fit^2, .Machine$double.eps * max(F2, 1))
    }
    bic[i] <- log(rss / np) + df * log(np) / np
  }
  best <- which.min(bic)
  list(lambda = lambdas[best], lambdas = lambdas, bic = bic, index = best)
}

# Rank-1 SFPCA in the constraint form of Allen & Weylandt (2019):
# alternate exact solves of
#   min_x 0.5 x'Sx - b'x + P(x; lambda),  S = I + alpha * Omega,
# (via C++ coordinate descent) with rescaling onto the S-norm ball.
# For convex P (l1) each half-step solves the constrained subproblem exactly,
# so the objective u'Xv - P(u) - P(v) is monotonically non-decreasing.
#' @noRd
sfpca_rank1 <- function(X = NULL,
                        lambda_u, lambda_v,
                        alpha_u, alpha_v,
                        Omega_u, Omega_v,
                        penalty_u, penalty_v,
                        max_iter, tol, verbose,
                        u_init = NULL, v_init = NULL, d_init = NULL,
                        ops = NULL, exact_inner = FALSE) {
  if (is.null(ops)) {
    if (is.null(X)) stop("Either X or ops must be supplied.")
    ops <- sfpca_make_ops(X)
  }
  n <- ops$n
  p <- ops$p
  # Initialize u and v
  if (is.null(u_init) || is.null(v_init)) {
    if (is.null(X)) stop("u_init and v_init are required when only ops is supplied.")
    svd_res <- svd1_deflated(X, verbose = verbose)
    u_init <- as.numeric(svd_res$u)
    v_init <- as.numeric(svd_res$v)
    d_init <- svd_res$d[1]
  }
  if (is.null(d_init)) {
    d_init <- sum(u_init * ops$mv(v_init))
  }

  S_u <- as_dgc(Matrix::Diagonal(n) + alpha_u * Omega_u)
  S_v <- as_dgc(Matrix::Diagonal(p) + alpha_v * Omega_v)

  # Start from the SVD initializer, normalized in the respective S-norms.
  u <- u_init / s_norm(u_init, S_u)
  v <- v_init / s_norm(v_init, S_v)
  # Warm starts for the subproblem solver live on the unnormalized scale
  # (the penalized solution scales with the current singular value).
  u_hat <- abs(d_init) * u_init
  v_hat <- abs(d_init) * v_init

  iter <- 0
  obj_old <- -Inf
  obj_trace <- numeric(0)
  kkt_u <- NA_real_
  kkt_v <- NA_real_
  degenerate <- FALSE
  conv_u <- TRUE
  conv_v <- TRUE
  # Inexact alternating minimization: early subproblem solves only need to be
  # as accurate as the outer progress (the target b moves anyway). The inner
  # KKT tolerance tracks the relative objective change and is forced down to
  # `inner_floor` before the outer loop may declare convergence, so the final
  # solves are always exact to tolerance. `exact_inner = TRUE` keeps every
  # solve at the floor (restores the monotone-objective guarantee for l1).
  inner_floor <- 1e-8
  inner_tol <- if (exact_inner) inner_floor else 1e-3
  repeat {
    iter <- iter + 1
    # --- u-update: solve of 0.5 u'S_u u - u'(Xv) + P(u; lambda_u)
    b_u <- ops$mv(v)
    sol_u <- sfpca_cd_solve(S_u, b_u, u_hat, lambda_u, penalty_u,
                            tol = inner_tol)
    u_hat <- sol_u$x
    kkt_u <- sol_u$kkt
    conv_u <- sol_u$converged
    if (all(u_hat == 0)) {
      degenerate <- TRUE
      break
    }
    u <- u_hat / s_norm(u_hat, S_u)

    # --- v-update: solve of 0.5 v'S_v v - v'(X'u) + P(v; lambda_v)
    b_v <- ops$tmv(u)
    sol_v <- sfpca_cd_solve(S_v, b_v, v_hat, lambda_v, penalty_v,
                            tol = inner_tol)
    v_hat <- sol_v$x
    kkt_v <- sol_v$kkt
    conv_v <- sol_v$converged
    if (all(v_hat == 0)) {
      degenerate <- TRUE
      break
    }
    v <- v_hat / s_norm(v_hat, S_v)

    # Objective of the constrained problem (monotone for convex penalties).
    # u'Xv = v'(X'u) reuses b_v: no extra matrix-vector product.
    obj <- sum(v * b_v) -
      penalty_value(u, penalty_u, lambda_u) -
      penalty_value(v, penalty_v, lambda_v)
    obj_trace[iter] <- obj
    if (verbose) cat("Iteration", iter, "Objective:", obj, "\n")

    # Check convergence
    stalled <- abs(obj - obj_old) < tol
    if (iter >= max_iter) break
    if (stalled) {
      if (inner_tol <= inner_floor) break
      # Progress stalled under loose inner solves: tighten to the floor so
      # the terminal iterates come from exact-to-tolerance subproblem solves.
      inner_tol <- inner_floor
    } else if (!exact_inner) {
      rel_change <- abs(obj - obj_old) / max(1, abs(obj))
      inner_tol <- max(inner_floor, min(1e-3, 0.1 * rel_change))
    }
    obj_old <- obj
  }

  if (!conv_u || !conv_v) {
    warning("Coordinate descent did not reach tolerance in the final ",
            "subproblem solves (KKT residuals: u = ", signif(kkt_u, 3),
            ", v = ", signif(kkt_v, 3), "). The smoothness penalty ",
            "S = I + alpha * Omega may be badly conditioned; consider a ",
            "smaller alpha or a rescaled Omega.")
  }

  if (degenerate) {
    # Penalty strong enough to zero out a factor: dead component.
    if (verbose) cat("Component is fully sparse (zero solution); returning d = 0.\n")
    return(list(u = rep(0, n), v = rep(0, p), d = 0,
                obj_trace = obj_trace, kkt_u = kkt_u, kkt_v = kkt_v,
                iters = iter))
  }

  # Final sign fix and Euclidean normalization
  if (sum(u) < 0) {
    u <- -u
    v <- -v
  }
  u <- u / sqrt(sum(u^2))
  v <- v / sqrt(sum(v^2))
  d <- sum(u * ops$mv(v))
  return(list(u = u, v = v, d = d,
              obj_trace = obj_trace, kkt_u = kkt_u, kkt_v = kkt_v,
              iters = iter))
}

# S-norm sqrt(x'Sx) with a floor to avoid division by zero
#' @noRd
s_norm <- function(x, S) {
  sqrt(max(as.numeric(Matrix::crossprod(x, S %*% x)), .Machine$double.eps))
}

# R wrapper around the C++ coordinate-descent subproblem solver.
# S must already be a dgCMatrix. The convergence tolerance is on the KKT
# residual, which lives on the gradient scale and grows with ||b||.
#' @noRd
sfpca_cd_solve <- function(S, b, x0, lambda, penalty,
                           max_sweeps = 1000L, tol = 1e-8, scad_a = 3.7) {
  code <- switch(penalty,
                 l1 = 0L,
                 scad = 1L,
                 stop("Unsupported penalty function: ", penalty))
  tol_abs <- tol * max(1, max(abs(b)))
  sfpca_cd_solve_cpp(S, as.numeric(b), as.numeric(x0), lambda, code,
                     scad_a, as.integer(max_sweeps), tol_abs)
}

# Penalty value P(x; lambda) (lambda included)
#' @noRd
penalty_value <- function(x, penalty, lambda, a = 3.7) {
  if (penalty == "l1") {
    lambda * sum(abs(x))
  } else if (penalty == "scad") {
    abs_x <- abs(x)
    penalty_values <- ifelse(
      abs_x <= lambda,
      lambda * abs_x,
      ifelse(
        abs_x <= a * lambda,
        (2 * a * lambda * abs_x - abs_x^2 - lambda^2) / (2 * (a - 1)),
        (a + 1) * lambda^2 / 2
      )
    )
    sum(penalty_values)
  } else {
    stop("Unsupported penalty function.")
  }
}
