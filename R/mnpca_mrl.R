#' Matrix-Normal PCA via Maximum Regularized Likelihood
#'
#' Fits a rank-\code{ncomp} matrix factorization under matrix-normal noise with
#' sparse row/column precision matrices:
#' \deqn{
#' Y = XW^\top + E,\quad E \sim \mathcal{MN}(0,\Omega,\Sigma)
#' }
#' using alternating block updates:
#' \enumerate{
#'   \item Alternating least squares updates for \code{X, W} with fixed
#'   precisions (\code{Theta_row = Omega^{-1}}, \code{Theta_col = Sigma^{-1}}).
#'   \item Graphical-lasso style precision updates for \code{Theta_row} and
#'   \code{Theta_col} using ADMM, warm starts, and optional block screening.
#' }
#' The precision updates are solved on a residual-free scatter surrogate
#' (rather than the exact residual \code{E = Y - XW^T}), so the two block
#' updates do not jointly minimize a single shared objective at each step;
#' as a result the reported \code{objective_path} is a convergence
#' diagnostic, not a monotone objective trace (see Details and the
#' \code{objective_path} return value).
#'
#' The covariance updates use low-rank correction identities and avoid explicit
#' construction of \code{E = Y - XW^T}.
#'
#' This function returns a plain S3 list of class \code{"mnpca_mrl"}, not a
#' \pkg{multivarious} \code{bi_projector}/\code{cross_projector}; the
#' \code{scores()}/\code{components()}/\code{reconstruct()} generics used
#' elsewhere in this package do not apply here. Use the returned
#' \code{$X}, \code{$W}, and \code{$fitted} fields directly.
#'
#' @param Y Numeric matrix (\code{n x p}).
#' @param ncomp Target rank (\code{r}).
#' @param lambda_row L1 penalty for row precision (\code{Theta_row}).
#' @param lambda_col L1 penalty for column precision (\code{Theta_col}).
#' @param max_outer Maximum number of outer BCD iterations.
#' @param max_inner Maximum ALS steps per outer iteration.
#' @param tol Relative tolerance used for ALS and objective convergence checks.
#' @param eps_ridge Ridge added to small \code{r x r} normal-equation systems.
#' @param jitter Small diagonal jitter used in covariance/precision updates.
#' @param center Logical; center columns of \code{Y} before fitting.
#' @param update_precisions Logical; if \code{FALSE}, keeps identity precisions
#'   and runs weighted ALS only.
#' @param warm_start Logical; warm start precision updates from previous iterate.
#' @param gl_maxit Maximum ADMM iterations per graphical-lasso subproblem.
#' @param gl_tol ADMM convergence tolerance for graphical-lasso subproblems.
#' @param gl_rho ADMM augmented Lagrangian parameter.
#' @param penalize_diagonal Logical; whether to penalize diagonal precision
#'   entries in L1 term. Default \code{FALSE}.
#' @param block_screen Logical; use thresholded connected components to solve
#'   precision subproblems blockwise.
#' @param scale_fix One of \code{"trace"} or \code{"none"}. \code{"trace"}
#'   normalizes each precision to mean diagonal 1 after updates.
#' @param sparsify_threshold Off-diagonal magnitude threshold used to set tiny
#'   precision entries to zero after each graphical-lasso solve.
#' @param as_sparse_precision Logical; store precisions as sparse matrices when
#'   many entries are zero.
#' @param verbose Logical; print iteration diagnostics.
#'
#' @return An object of class \code{"mnpca_mrl"} with components:
#' \describe{
#'   \item{X, W}{Estimated low-rank factors (\code{n x r}, \code{p x r}).}
#'   \item{Theta_row, Theta_col}{Estimated row/column precision matrices.}
#'   \item{fitted}{Reconstructed matrix on input scale.}
#'   \item{fitted_centered}{Reconstructed centered matrix used in optimization.}
#'   \item{residual_centered}{Centered residual matrix.}
#'   \item{objective_path}{Objective value evaluated after each outer
#'     iteration's block updates, before any trace rescaling of the
#'     precisions. Because the factor updates (ALS) and precision updates
#'     (graphical lasso on a residual-free scatter surrogate) target
#'     different surrogates rather than a single shared objective, this
#'     path is a convergence heuristic and is \strong{not guaranteed to be
#'     monotone}.}
#'   \item{iterations}{Number of outer iterations used.}
#'   \item{converged}{Logical; \code{TRUE} if the relative change in
#'     \code{objective_path} fell below \code{tol} before \code{max_outer}
#'     was reached. This is a convergence heuristic based on relative
#'     objective change, not a guarantee of a local optimum.}
#'   \item{center}{Column centering vector (or \code{NULL}).}
#'   \item{call}{Matched call.}
#' }
#'
#' @details
#' This implementation follows the maximum regularized likelihood (MRL)
#' formulation of MN-PCA, combining low-rank factor updates with sparse
#' precision estimation in row and column spaces, via alternating surrogate
#' updates rather than joint minimization of a single objective at every
#' step.
#'
#' The main optimization target is:
#' \deqn{
#' \frac12\mathrm{tr}\left(\Theta_c (Y-XW^\top)^\top \Theta_r (Y-XW^\top)\right)
#' -\frac{p}{2}\log|\Theta_r|-\frac{n}{2}\log|\Theta_c|
#' + n\lambda_r\|\Theta_r\|_1 + p\lambda_c\|\Theta_c\|_1
#' }
#' with \eqn{\Theta_r \succ 0, \Theta_c \succ 0}. When
#' \code{update_precisions = FALSE}, the method reduces to weighted low-rank
#' approximation with fixed identity precisions.
#'
#' Because the ALS factor updates and the graphical-lasso precision updates
#' are solved against different surrogates of this target (the precision
#' step uses a residual-free scatter matrix rather than the exact residual),
#' the outer alternation is best understood as alternating surrogate updates
#' with a relative-objective-change convergence heuristic, not classical
#' block coordinate descent on a single monotone objective.
#'
#' @references
#' Zhang, C., Gai, K., & Zhang, S. (2024).
#' \emph{Matrix normal PCA for interpretable dimension reduction and graphical noise modeling}.
#' Pattern Recognition, 154, 110591.
#' \doi{10.1016/j.patcog.2024.110591}
#'
#' Friedman, J., Hastie, T., & Tibshirani, R. (2008).
#' \emph{Sparse inverse covariance estimation with the graphical lasso}.
#' Biostatistics, 9(3), 432-441.
#'
#' @examples
#' set.seed(123)
#' n <- 20; p <- 12; r <- 3
#' Y <- matrix(rnorm(n * p), n, p)
#' fit <- mnpca_mrl(
#'   Y,
#'   ncomp = r,
#'   lambda_row = 0.08,
#'   lambda_col = 0.08,
#'   max_outer = 6,
#'   max_inner = 6,
#'   verbose = FALSE
#' )
#' dim(fit$X)
#' dim(fit$W)
#' length(fit$objective_path)
#'
#' @export
mnpca_mrl <- function(Y,
                      ncomp = min(dim(Y)),
                      lambda_row = 0.05,
                      lambda_col = 0.05,
                      max_outer = 25,
                      max_inner = 5,
                      tol = 1e-4,
                      eps_ridge = 1e-8,
                      jitter = 1e-6,
                      center = TRUE,
                      update_precisions = TRUE,
                      warm_start = TRUE,
                      gl_maxit = 200,
                      gl_tol = 1e-4,
                      gl_rho = 1,
                      penalize_diagonal = FALSE,
                      block_screen = TRUE,
                      scale_fix = c("trace", "none"),
                      sparsify_threshold = 1e-8,
                      as_sparse_precision = TRUE,
                      verbose = FALSE) {
  call <- match.call()
  scale_fix <- match.arg(scale_fix)

  if (inherits(Y, "Matrix")) {
    Y <- as.matrix(Y)
  } else {
    Y <- as.matrix(Y)
  }
  if (!is.numeric(Y)) {
    stop("Y must be numeric.")
  }
  if (any(!is.finite(Y))) {
    stop("Y contains non-finite values.")
  }

  n <- nrow(Y)
  p <- ncol(Y)
  if (n == 0L || p == 0L) {
    stop("Y must have positive dimensions.")
  }

  ncomp <- as.integer(ncomp)
  if (!is.finite(ncomp) || ncomp < 1L || ncomp > min(n, p)) {
    stop("ncomp must be an integer in [1, min(nrow(Y), ncol(Y))].")
  }
  if (max_outer < 1L || max_inner < 1L) {
    stop("max_outer and max_inner must be >= 1.")
  }
  if (lambda_row < 0 || lambda_col < 0) {
    stop("lambda_row and lambda_col must be non-negative.")
  }
  if (eps_ridge <= 0 || jitter <= 0) {
    stop("eps_ridge and jitter must be positive.")
  }

  center_vec <- NULL
  Y_work <- Y
  if (isTRUE(center)) {
    center_vec <- colMeans(Y_work)
    Y_work <- sweep(Y_work, 2L, center_vec, FUN = "-")
  }

  init <- .mnpca_init_factors(Y_work, ncomp)
  X <- init$X
  W <- init$W

  Theta_row <- Matrix::Diagonal(n)
  Theta_col <- Matrix::Diagonal(p)

  objective_path <- numeric(0)
  prev_obj <- Inf
  converged <- FALSE
  outer_used <- 0L

  for (it in seq_len(max_outer)) {
    outer_used <- it

    for (k in seq_len(max_inner)) {
      A <- .mnpca_spmm(Theta_col, W)
      G <- Matrix::crossprod(W, A) + eps_ridge * diag(ncomp)
      Ginv <- .mnpca_safe_solve_spd(G, eps = eps_ridge)
      P <- Y_work %*% A
      X_new <- P %*% Ginv

      B <- .mnpca_spmm(Theta_row, X_new)
      H <- Matrix::crossprod(X_new, B) + eps_ridge * diag(ncomp)
      Hinv <- .mnpca_safe_solve_spd(H, eps = eps_ridge)
      Q <- Matrix::t(Y_work) %*% B
      W_new <- Q %*% Hinv

      dx <- .mnpca_relative_change(X_new, X)
      dw <- .mnpca_relative_change(W_new, W)
      X <- X_new
      W <- W_new
      if (dx < tol && dw < tol) {
        break
      }
    }

    if (isTRUE(update_precisions)) {
      scat <- .mnpca_scatter_no_residual(
        Y = Y_work,
        X = X,
        W = W,
        Theta_row = Theta_row,
        Theta_col = Theta_col,
        jitter = jitter
      )

      row_init <- if (isTRUE(warm_start)) Theta_row else NULL
      col_init <- if (isTRUE(warm_start)) Theta_col else NULL

      row_fit <- .mnpca_glasso_admm(
        S = scat$S1,
        lambda = lambda_row,
        Theta_init = row_init,
        rho = gl_rho,
        maxit = gl_maxit,
        tol = gl_tol,
        penalize_diagonal = penalize_diagonal,
        jitter = jitter,
        block_screen = block_screen
      )
      col_fit <- .mnpca_glasso_admm(
        S = scat$S2,
        lambda = lambda_col,
        Theta_init = col_init,
        rho = gl_rho,
        maxit = gl_maxit,
        tol = gl_tol,
        penalize_diagonal = penalize_diagonal,
        jitter = jitter,
        block_screen = block_screen
      )

      Theta_row <- .mnpca_postprocess_precision(
        row_fit$Theta,
        jitter = jitter,
        threshold = sparsify_threshold,
        as_sparse = as_sparse_precision
      )
      Theta_col <- .mnpca_postprocess_precision(
        col_fit$Theta,
        jitter = jitter,
        threshold = sparsify_threshold,
        as_sparse = as_sparse_precision
      )

    }

    # Evaluate the objective on the iterates the block updates actually
    # produced, BEFORE any identifiability rescaling: the matrix-normal
    # likelihood is invariant only under *joint* reciprocal rescaling
    # (c*Theta_row, Theta_col/c), so two independent trace normalizations
    # change the objective and would make the reported path non-monotone
    # (and the convergence test unreliable).
    obj <- .mnpca_objective_value(
      Y = Y_work,
      X = X,
      W = W,
      Theta_row = Theta_row,
      Theta_col = Theta_col,
      lambda_row = lambda_row,
      lambda_col = lambda_col,
      penalize_diagonal = penalize_diagonal
    )
    objective_path <- c(objective_path, obj)

    if (isTRUE(verbose)) {
      message(sprintf("[mnpca_mrl] iter=%d obj=%.6f", it, obj))
    }

    if (it > 1L) {
      rel_obj <- abs(prev_obj - obj) / max(1, abs(prev_obj))
      if (rel_obj < tol) {
        converged <- TRUE
        break
      }
    }
    prev_obj <- obj

    # Identifiability projection between iterations (prevents scale drift in
    # the glasso subproblems). Applied after the objective/convergence test.
    if (isTRUE(update_precisions) && scale_fix == "trace") {
      tr_row <- mean(Matrix::diag(Theta_row))
      tr_col <- mean(Matrix::diag(Theta_col))
      if (is.finite(tr_row) && tr_row > 0) {
        Theta_row <- Theta_row / tr_row
      }
      if (is.finite(tr_col) && tr_col > 0) {
        Theta_col <- Theta_col / tr_col
      }
    }
  }

  # Returned precisions honor scale_fix regardless of which iteration broke
  # the loop.
  if (isTRUE(update_precisions) && scale_fix == "trace") {
    tr_row <- mean(Matrix::diag(Theta_row))
    tr_col <- mean(Matrix::diag(Theta_col))
    if (is.finite(tr_row) && tr_row > 0) {
      Theta_row <- Theta_row / tr_row
    }
    if (is.finite(tr_col) && tr_col > 0) {
      Theta_col <- Theta_col / tr_col
    }
  }

  # Canonicalize factors so that X^T Theta_row X ~= I.
  C <- Matrix::crossprod(X, .mnpca_spmm(Theta_row, X))
  C <- .mnpca_symmetrize(as.matrix(C)) + eps_ridge * diag(ncomp)
  R <- chol(C)
  Rinv <- backsolve(R, diag(ncomp))
  X <- X %*% Rinv
  W <- W %*% Matrix::t(R)

  fitted_centered <- X %*% Matrix::t(W)
  residual_centered <- Y_work - fitted_centered
  fitted <- fitted_centered
  if (!is.null(center_vec)) {
    fitted <- sweep(fitted, 2L, center_vec, FUN = "+")
  }

  ret <- list(
    X = X,
    W = W,
    Theta_row = Theta_row,
    Theta_col = Theta_col,
    fitted = fitted,
    fitted_centered = fitted_centered,
    residual_centered = residual_centered,
    objective_path = objective_path,
    iterations = outer_used,
    converged = converged,
    center = center_vec,
    call = call
  )
  class(ret) <- "mnpca_mrl"
  ret
}

.mnpca_scatter_no_residual <- function(Y, X, W, Theta_row, Theta_col, jitter = 0) {
  n <- nrow(Y)
  p <- ncol(Y)

  A <- .mnpca_spmm(Theta_col, W)
  G <- Matrix::crossprod(W, A)
  P <- Y %*% A

  B <- .mnpca_spmm(Theta_row, X)
  H <- Matrix::crossprod(X, B)
  Q <- Matrix::t(Y) %*% B

  TY <- .mnpca_spmm(Theta_col, Matrix::t(Y))
  M1 <- Y %*% TY
  R1 <- M1 - (P %*% Matrix::t(X)) - (X %*% Matrix::t(P)) + (X %*% G %*% Matrix::t(X))
  R1 <- .mnpca_symmetrize(R1)
  S1 <- R1 / p
  if (jitter > 0) {
    diag(S1) <- diag(S1) + jitter
  }

  UY <- .mnpca_spmm(Theta_row, Y)
  M2 <- Matrix::t(Y) %*% UY
  R2 <- M2 - (Q %*% Matrix::t(W)) - (W %*% Matrix::t(Q)) + (W %*% H %*% Matrix::t(W))
  R2 <- .mnpca_symmetrize(R2)
  S2 <- R2 / n
  if (jitter > 0) {
    diag(S2) <- diag(S2) + jitter
  }

  list(
    S1 = S1,
    S2 = S2,
    R1 = R1,
    R2 = R2,
    A = A,
    B = B,
    G = G,
    H = H,
    P = P,
    Q = Q
  )
}

.mnpca_glasso_admm <- function(S, lambda, Theta_init = NULL, rho = 1,
                               maxit = 200, tol = 1e-4,
                               penalize_diagonal = FALSE,
                               jitter = 1e-8,
                               block_screen = TRUE) {
  S <- .mnpca_symmetrize(as.matrix(S))
  p <- nrow(S)
  if (p != ncol(S)) {
    stop("S must be square.")
  }
  if (p == 1L) {
    sval <- max(S[1, 1] + jitter, jitter)
    denom <- if (penalize_diagonal) sval + lambda else sval
    return(list(Theta = matrix(1 / denom, 1, 1), converged = TRUE, iterations = 1L))
  }
  diag(S) <- diag(S) + jitter

  if (lambda <= 0) {
    Theta <- .mnpca_symmetrize(.mnpca_safe_solve_spd(S, eps = jitter))
    return(list(Theta = Theta, converged = TRUE, iterations = 1L))
  }

  if (!isTRUE(block_screen)) {
    return(.mnpca_glasso_admm_core(
      S = S,
      lambda = lambda,
      Theta_init = Theta_init,
      rho = rho,
      maxit = maxit,
      tol = tol,
      penalize_diagonal = penalize_diagonal,
      jitter = jitter
    ))
  }

  comps <- .mnpca_connected_components(S, lambda)
  if (length(comps) <= 1L) {
    return(.mnpca_glasso_admm_core(
      S = S,
      lambda = lambda,
      Theta_init = Theta_init,
      rho = rho,
      maxit = maxit,
      tol = tol,
      penalize_diagonal = penalize_diagonal,
      jitter = jitter
    ))
  }

  Theta <- matrix(0, p, p)
  converged_all <- TRUE
  max_iter <- 1L
  for (idx in comps) {
    init_block <- NULL
    if (!is.null(Theta_init)) {
      init_block <- as.matrix(Theta_init)[idx, idx, drop = FALSE]
    }
    block_fit <- .mnpca_glasso_admm_core(
      S = S[idx, idx, drop = FALSE],
      lambda = lambda,
      Theta_init = init_block,
      rho = rho,
      maxit = maxit,
      tol = tol,
      penalize_diagonal = penalize_diagonal,
      jitter = jitter
    )
    Theta[idx, idx] <- block_fit$Theta
    converged_all <- converged_all && isTRUE(block_fit$converged)
    max_iter <- max(max_iter, block_fit$iterations)
  }

  Theta <- .mnpca_symmetrize(Theta)
  diag(Theta) <- pmax(diag(Theta), jitter)
  list(Theta = Theta, converged = converged_all, iterations = max_iter)
}

.mnpca_glasso_admm_core <- function(S, lambda, Theta_init = NULL, rho = 1,
                                    maxit = 200, tol = 1e-4,
                                    penalize_diagonal = FALSE,
                                    jitter = 1e-8) {
  p <- nrow(S)
  S <- .mnpca_symmetrize(as.matrix(S))
  diag(S) <- diag(S) + jitter

  if (is.null(Theta_init)) {
    Theta <- diag(1 / pmax(diag(S), jitter), p)
  } else {
    Theta <- .mnpca_symmetrize(as.matrix(Theta_init))
    diag(Theta) <- pmax(diag(Theta), jitter)
  }
  Z <- Theta
  U <- matrix(0, p, p)

  r_norm <- Inf
  s_norm <- Inf
  eps_pri <- Inf
  eps_dual <- Inf
  sqrt_pp <- sqrt(p * p)
  it <- 0L

  for (k in seq_len(maxit)) {
    it <- k
    eig <- eigen(rho * (Z - U) - S, symmetric = TRUE)
    d <- eig$values
    Q <- eig$vectors
    x <- (d + sqrt(d * d + 4 * rho)) / (2 * rho)
    Theta <- Q %*% (x * Matrix::t(Q))
    Theta <- .mnpca_symmetrize(Theta)

    Z_old <- Z
    V <- Theta + U
    if (lambda > 0) {
      Z <- .mnpca_soft_threshold(V, lambda / rho)
    } else {
      Z <- V
    }
    if (!isTRUE(penalize_diagonal)) {
      diag(Z) <- diag(V)
    }
    Z <- .mnpca_symmetrize(Z)

    U <- U + Theta - Z

    r_norm <- .mnpca_fro_norm(Theta - Z)
    s_norm <- rho * .mnpca_fro_norm(Z - Z_old)
    eps_pri <- sqrt_pp * tol + tol * max(.mnpca_fro_norm(Theta), .mnpca_fro_norm(Z))
    eps_dual <- sqrt_pp * tol + tol * rho * .mnpca_fro_norm(U)
    if (r_norm <= eps_pri && s_norm <= eps_dual) {
      break
    }
  }

  Theta <- .mnpca_symmetrize(Theta)
  diag(Theta) <- pmax(diag(Theta), jitter)
  Theta <- as.matrix(ensure_spd(Theta, tol = jitter))
  list(
    Theta = Theta,
    converged = (r_norm <= eps_pri && s_norm <= eps_dual),
    iterations = it
  )
}

.mnpca_postprocess_precision <- function(Theta, jitter = 1e-8, threshold = 0,
                                         as_sparse = TRUE) {
  Theta <- .mnpca_symmetrize(as.matrix(Theta))
  if (threshold > 0) {
    offdiag <- row(Theta) != col(Theta)
    Theta[offdiag & abs(Theta) < threshold] <- 0
  }
  diag(Theta) <- pmax(diag(Theta), jitter)
  Theta <- as.matrix(ensure_spd(Theta, tol = jitter))

  if (isTRUE(as_sparse)) {
    zero_rate <- mean(Theta == 0)
    if (is.finite(zero_rate) && zero_rate > 0.4) {
      Theta_sp <- Matrix::Matrix(Theta, sparse = TRUE)
      if (methods::is(Theta_sp, "symmetricMatrix")) {
        Theta_sp <- methods::as(Theta_sp, "generalMatrix")
      }
      Theta_sp <- methods::as(Theta_sp, "CsparseMatrix")
      return(Theta_sp)
    }
  }
  Theta
}

.mnpca_objective_value <- function(Y, X, W, Theta_row, Theta_col,
                                   lambda_row, lambda_col,
                                   penalize_diagonal = FALSE) {
  n <- nrow(Y)
  p <- ncol(Y)
  scat <- .mnpca_scatter_no_residual(
    Y = Y,
    X = X,
    W = W,
    Theta_row = Theta_row,
    Theta_col = Theta_col,
    jitter = 0
  )

  Theta_row_d <- as.matrix(Theta_row)
  Theta_col_d <- as.matrix(Theta_col)
  quad <- sum(Theta_row_d * scat$R1)

  logdet_row <- tryCatch(
    as.numeric(Matrix::determinant(Matrix::forceSymmetric(Matrix::Matrix(Theta_row_d, sparse = FALSE),
                                                          uplo = "U"),
                                   logarithm = TRUE)$modulus),
    error = function(e) NA_real_
  )
  logdet_col <- tryCatch(
    as.numeric(Matrix::determinant(Matrix::forceSymmetric(Matrix::Matrix(Theta_col_d, sparse = FALSE),
                                                          uplo = "U"),
                                   logarithm = TRUE)$modulus),
    error = function(e) NA_real_
  )

  if (!is.finite(logdet_row) || !is.finite(logdet_col) || !is.finite(quad)) {
    return(Inf)
  }

  pen <- n * lambda_row * .mnpca_l1_penalty(Theta_row_d, penalize_diagonal = penalize_diagonal) +
    p * lambda_col * .mnpca_l1_penalty(Theta_col_d, penalize_diagonal = penalize_diagonal)

  0.5 * quad - 0.5 * p * logdet_row - 0.5 * n * logdet_col + pen
}

.mnpca_init_factors <- function(Y, ncomp) {
  n <- nrow(Y)
  p <- ncol(Y)
  sv <- tryCatch({
    if (min(n, p) > ncomp + 1L) {
      .top_svd(Y, ncomp, nu = ncomp, nv = ncomp)
    } else {
      base::svd(Y, nu = ncomp, nv = ncomp)
    }
  }, error = function(e) {
    base::svd(Y, nu = ncomp, nv = ncomp)
  })

  U <- as.matrix(sv$u[, seq_len(ncomp), drop = FALSE])
  V <- as.matrix(sv$v[, seq_len(ncomp), drop = FALSE])
  d <- as.numeric(sv$d[seq_len(ncomp)])
  X <- U %*% diag(d, nrow = ncomp, ncol = ncomp)
  list(X = X, W = V)
}

.mnpca_connected_components <- function(S, lambda) {
  p <- nrow(S)
  if (p == 1L || lambda <= 0) {
    return(list(seq_len(p)))
  }
  A <- abs(S) > lambda
  diag(A) <- FALSE

  seen <- rep(FALSE, p)
  comps <- list()
  n_comp <- 0L
  for (i in seq_len(p)) {
    if (seen[i]) next
    n_comp <- n_comp + 1L
    stack <- i
    seen[i] <- TRUE
    comp <- integer(0)
    while (length(stack) > 0L) {
      v <- stack[[length(stack)]]
      stack <- stack[-length(stack)]
      comp <- c(comp, v)
      nbrs <- which(A[v, ] & !seen)
      if (length(nbrs) > 0L) {
        seen[nbrs] <- TRUE
        stack <- c(stack, nbrs)
      }
    }
    comps[[n_comp]] <- sort(comp)
  }
  comps
}

.mnpca_spmm <- function(Theta, X) {
  as.matrix(Theta %*% X)
}

.mnpca_safe_solve_spd <- function(G, eps = 1e-8) {
  G <- .mnpca_symmetrize(as.matrix(G))
  bump <- eps
  for (i in 0:7) {
    G_try <- G
    diag(G_try) <- diag(G_try) + bump
    out <- tryCatch({
      R <- chol(G_try)
      chol2inv(R)
    }, error = function(e) NULL)
    if (!is.null(out)) {
      return(out)
    }
    bump <- bump * 10
  }
  stop("Failed SPD solve in .mnpca_safe_solve_spd.")
}

.mnpca_symmetrize <- function(M) {
  0.5 * (M + Matrix::t(M))
}

.mnpca_soft_threshold <- function(M, tval) {
  sign(M) * pmax(abs(M) - tval, 0)
}

.mnpca_fro_norm <- function(M) {
  sqrt(sum(M * M))
}

.mnpca_relative_change <- function(new, old) {
  .mnpca_fro_norm(new - old) / max(1, .mnpca_fro_norm(old))
}

.mnpca_l1_penalty <- function(Theta, penalize_diagonal = FALSE) {
  if (isTRUE(penalize_diagonal)) {
    return(sum(abs(Theta)))
  }
  off <- Theta
  diag(off) <- 0
  sum(abs(off))
}
