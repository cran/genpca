#' Generalized PLS-SVD via Implicit Operator (memory-safe)
#'
#' Compute the top-k singular triplets of \eqn{S = Xe' Ye} without
#' materializing the whitened matrices \eqn{Xe = Mx^{1/2} X Wx^{1/2}},
#' \eqn{Ye = My^{1/2} Y Wy^{1/2}} when doing so would densify sparse data.
#' When the whitening is sparsity-preserving (identity/diagonal metrics) or
#' the data are dense, the whitened blocks are precomputed once so each
#' matrix-vector product in the iterative SVD costs two multiplies.
#'
#' Naming map: this function names its metrics `XLW`/`YLW` (left/row
#' weights) and `XRW`/`YRW` (right/column weights); these correspond to
#' `Mx`/`My` (row metrics) and `Ax`/`Ay` (column metrics) in `genpca()`'s and
#' `genpls()`'s M/A convention.
#'
#' @param X n x I matrix (numeric or Matrix)
#' @param Y n x J matrix (numeric or Matrix)
#' @param XLW Row metric for X (M_X): NULL/identity, numeric length-n, diagonalMatrix, or PSD Matrix
#' @param YLW Row metric for Y (M_Y)
#' @param XRW Column metric for X (W_X)
#' @param YRW Column metric for Y (W_Y)
#' @param k  Number of components. If `k` exceeds
#'   `min(ncol(X), ncol(Y))`, a warning is issued and `k` is silently
#'   truncated to that maximum.
#' @param center,scale Logical; pre-center/scale columns of X, Y before metrics
#' @param svd_backend One of "eigencore" (default) or "irlba"; "RSpectra" is
#'   accepted as a deprecated alias of "eigencore". Ignored
#'   whenever both `ncol(X) <= 64` and `ncol(Y) <= 64`, in which case `S` is
#'   materialized densely and solved with `base::svd()`.
#' @param svd_opts List of options for the backend: `tol` (both backends) and
#'   `maxitr` (irlba only; the eigencore partial SVD has no iteration cap).
#'   An incomplete eigencore solve raises `genpca_solver_nonconvergence`;
#'   try a less stringent `tol` if the requested accuracy cannot be reached.
#' @param constraints_remedy What to do with a metric that is not positive
#'   semi-definite: `"error"` (default), `"ridge"`, `"clip"` or `"identity"`;
#'   repairs emit a `genpca_metric_repaired` warning. See [genpca()].
#' @return A list with elements:
#'   \describe{
#'     \item{d}{Length-`k` numeric vector of singular values of
#'       \eqn{S = Xe' Ye}.}
#'     \item{u}{`I x k` matrix; left singular vectors of `S` (orthonormal in
#'       the Euclidean metric).}
#'     \item{v}{`J x k` matrix; right singular vectors of `S` (orthonormal in
#'       the Euclidean metric).}
#'     \item{p}{`I x k` matrix of generalized X-weights,
#'       \eqn{p = W_X^{-1/2} u}.}
#'     \item{q}{`J x k` matrix of generalized Y-weights,
#'       \eqn{q = W_Y^{-1/2} v}.}
#'     \item{fi}{`I x k` matrix of X-variable scores, \eqn{F_i = W_X p D}
#'       (columns of `p` rescaled by the singular values).}
#'     \item{fj}{`J x k` matrix of Y-variable scores, \eqn{F_j = W_Y q D}.}
#'     \item{lx}{`N x k` matrix of X row latent variables,
#'       \eqn{L_x = M_X^{1/2} X W_X p}.}
#'     \item{ly}{`N x k` matrix of Y row latent variables,
#'       \eqn{L_y = M_Y^{1/2} Y W_Y q}.}
#'     \item{k}{Integer; number of components actually returned (may be
#'       less than the requested `k` if it exceeded `min(I, J)`).}
#'     \item{dims}{A list `list(N, I, J)` with the row count `N` and column
#'       counts `I = ncol(X)`, `J = ncol(Y)`.}
#'     \item{center}{A list `list(X, Y)` of the length-`I` / length-`J`
#'       column means subtracted from `X`/`Y` (all zero when
#'       `center = FALSE`).}
#'     \item{scale}{A list `list(X, Y)` of the length-`I` / length-`J`
#'       column scale factors divided out of `X`/`Y` (all one when
#'       `scale = FALSE`).}
#'   }
#' @references
#' Beaton, D. (2020). Generalized eigen, singular value, and partial least
#' squares decompositions: The GSVD package. arXiv:2010.14734.
#'
#' Abdi, H. (2007). Partial least square regression PLS-Regression. In
#' N. Salkind (Ed.), *Encyclopedia of Measurement and Statistics*.
#' Thousand Oaks, CA: Sage.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(40 * 6), 40, 6)
#' Y <- matrix(rnorm(40 * 4), 40, 4)
#' op <- gplssvd_op(X, Y, k = 2, center = TRUE)
#' round(op$d, 3)
#' @export
gplssvd_op <- function(X, Y,
                       XLW = NULL, YLW = NULL,
                       XRW = NULL, YRW = NULL,
                       k = 2, center = FALSE, scale = FALSE,
                       svd_backend = c("eigencore", "irlba", "RSpectra"),
                       svd_opts = list(tol = 1e-7, maxitr = 1000),
                       constraints_remedy = c("error", "ridge", "clip", "identity")) {

  svd_backend <- match.arg(svd_backend)
  constraints_remedy <- match.arg(constraints_remedy)
  if (svd_backend == "RSpectra") svd_backend <- "eigencore"
  if (!is.numeric(k) || length(k) != 1 || k < 1) {
    stop("k must be a single positive integer >= 1")
  }
  k <- as.integer(k)

  N <- nrow(X)
  I <- ncol(X)
  J <- ncol(Y)
  stopifnot(nrow(Y) == N)

  # Validate k against matrix dimensions
  if (k > min(I, J)) {
    warning("k (", k, ") exceeds min(ncol(X), ncol(Y)) = ", min(I, J),
            "; will return at most ", min(I, J), " components")
    k <- min(I, J)
  }

  # Column center/scale prior to constraints; keeps base matrices base and
  # Matrix objects Matrix (centering a sparse matrix necessarily densifies).
  cs <- function(A, do_center, do_scale) {
    cen <- rep(0, ncol(A))
    scl <- rep(1, ncol(A))
    if (isTRUE(do_center)) {
      cen <- as.numeric(Matrix::colMeans(A))
      if (inherits(A, "Matrix")) {
        A <- A - Matrix::Matrix(rep(1, nrow(A)), ncol = 1) %*%
          Matrix::Matrix(cen, nrow = 1)
      } else {
        A <- A - rep(cen, each = nrow(A))
      }
    }
    if (isTRUE(do_scale)) {
      s <- sqrt(Matrix::colSums(A^2) / pmax(nrow(A) - 1, 1))
      s[s == 0] <- 1
      scl <- as.numeric(s)
      if (inherits(A, "Matrix")) {
        A <- A %*% Matrix::Diagonal(x = 1 / scl)
      } else {
        A <- A * rep(1 / scl, each = nrow(A))
      }
    }
    list(A = A, center = cen, scale = scl)
  }

  Xcs <- cs(X, center, scale)
  Ycs <- cs(Y, center, scale)
  X <- Xcs$A
  Y <- Ycs$A

  # Metric operators (shared helper); reuse when both blocks share a metric
  MX <- .metric_operators(XLW, N, remedy = constraints_remedy, name = "XLW")
  MY <- if (identical(XLW, YLW)) MX else .metric_operators(YLW, N, remedy = constraints_remedy, name = "YLW")
  WX <- .metric_operators(XRW, I, remedy = constraints_remedy, name = "XRW")
  WY <- if (identical(XRW, YRW) && (is.null(XRW) || J == I)) {
    WX
  } else {
    .metric_operators(YRW, J, remedy = constraints_remedy, name = "YRW")
  }

  # Linear operators for S = t(Xe) %*% Ye (shared builder)
  opc <- .build_pls_operator(X, Y, MX, MY, WX, WY)

  # Small-dense fallback for stability on toy sizes
  use_dense <- (I <= 64 && J <= 64)
  if (!use_dense) {
    if (svd_backend == "irlba" && !requireNamespace("irlba", quietly = TRUE)) {
      stop("irlba package required for svd_backend='irlba'.")
    }
  }

  if (use_dense) {
    # Build S explicitly in the same algebra as the implicit operator:
    # S = WX^{1/2} X^T MX^{1/2} MY^{1/2} Y WY^{1/2}
    S <- if (isTRUE(opc$materialized)) {
      Matrix::crossprod(opc$Xe, opc$Ye)
    } else {
      T3 <- Matrix::crossprod(X, MX$mult_sqrt(MY$mult_sqrt(Y)))
      WY$mult_sqrt_right(WX$mult_sqrt(T3))
    }
    svdS <- svd(as.matrix(S))
    u <- svdS$u[, seq_len(k), drop = FALSE]
    v <- svdS$v[, seq_len(k), drop = FALSE]
    d <- svdS$d[seq_len(k)]
  } else if (svd_backend == "eigencore") {
    sv <- .top_svd(opc$S_mv, k, nu = k, nv = k,
                   tol = if (!is.null(svd_opts$tol)) svd_opts$tol else 1e-8,
                   adjoint = opc$ST_mv, dim = c(I, J))
    u <- sv$u
    v <- sv$v
    d <- sv$d
  } else {
    # irlba operator path via mult() callback. The dummy A supplies dims only,
    # so use an empty sparse matrix rather than allocating I x J dense zeros.
    A0 <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                               dims = c(I, J), x = numeric(0))
    mult_fun <- function(x, y) {
      # Handle both mult(A, v) and mult(v, A) calling styles
      if (is.matrix(x) || inherits(x, "Matrix")) {
        # x is A, y is vector(s): return A %*% y
        opc$S_mv(y)
      } else {
        # x is vector(s), y is A: return t(A) %*% x
        opc$ST_mv(x)
      }
    }
    sv <- irlba::irlba(A0,
                       nv = k, nu = k,
                       work = max(3 * k, k + 1),
                       tol = if (!is.null(svd_opts$tol)) svd_opts$tol else 1e-7,
                       maxit = if (!is.null(svd_opts$maxitr)) svd_opts$maxitr else 1000,
                       mult = mult_fun,
                       fastpath = FALSE)
    u <- sv$u
    v <- sv$v
    d <- sv$d
  }

  u <- as.matrix(u)
  v <- as.matrix(v)

  # Generalized singular vectors & component scores
  p <- as.matrix(WX$mult_invsqrt(u))
  q <- as.matrix(WY$mult_invsqrt(v))
  # Fi = WX p D, Fj = WY q D: `rep(d, each = nrow)` scales column j by d[j]
  Fi <- as.matrix(WX$mult(p))
  Fi <- Fi * rep(d, each = nrow(Fi))
  Fj <- as.matrix(WY$mult(q))
  Fj <- Fj * rep(d, each = nrow(Fj))
  # Latent variables Lx = MX^{1/2} X WX p = Xe u (identical for symmetric
  # PSD square roots, since WX WX^{-1/2} = WX^{1/2} on the range of WX)
  if (isTRUE(opc$materialized)) {
    Lx <- as.matrix(opc$Xe %*% u)
    Ly <- as.matrix(opc$Ye %*% v)
  } else {
    Lx <- as.matrix(MX$mult_sqrt(X %*% WX$mult(p)))
    Ly <- as.matrix(MY$mult_sqrt(Y %*% WY$mult(q)))
  }

  list(
    d  = as.numeric(d),
    u  = u,
    v  = v,
    p  = p,
    q  = q,
    fi = Fi,
    fj = Fj,
    lx = Lx,
    ly = Ly,

    k = k,
    dims = list(N = N, I = I, J = J),
    center = list(X = Xcs$center, Y = Ycs$center),
    scale  = list(X = Xcs$scale,  Y = Ycs$scale)
  )
}
