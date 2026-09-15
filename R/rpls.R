#' @noRd
fit_rpls <- function(X, Y,
                     K          = 2,
                     lambda     = 0.1,
                     penalty    = c("l1", "ridge"),
                     Q          = NULL,
                     nonneg     = FALSE,
                     tol        = 1e-6,
                     maxiter    = 200,
                     verbose    = FALSE) {
  penalty <- match.arg(penalty)

  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Y)
  if (n != nrow(Y)) {
    stop("X and Y must have the same number of rows")
  }

  if (!is.numeric(lambda) || any(is.na(lambda))) {
    stop("lambda must be numeric and non-missing")
  }
  if (any(lambda < 0)) {
    stop("lambda must be non-negative")
  }

  # Make lambda into a length-K vector
  if (length(lambda) == 1L) {
    lambda.vec <- rep(lambda, K)
  } else if (length(lambda) < K) {
    stop("lambda must be either a single numeric or length >= K.")
  } else {
    lambda.vec <- lambda
  }

  # If Q is null => identity, else must be p x p
  Id <- diag(p) # Keep identity for later use if needed
  if (is.null(Q)) {
    Q <- Id
    use.gpls <- FALSE
  } else {
    # Accept sparse or dense, etc.
    # Check dimension & basic positivity
    if (!all(dim(Q) == c(p, p))) {
      stop("Q must be p x p.")
    }
    # TODO: Check positive semi-definiteness? Might be slow.
    use.gpls <- TRUE
  }

  # Cross-product matrix
  M <- crossprod(X, Y)

  # Allocate
  V <- matrix(0, p, K) # X-loadings
  U <- matrix(0, q, K) # Y-loadings
  Z <- matrix(0, n, K) # X-scores (factors)

  # For SIMPLS deflation
  Rk <- matrix(0, p, 0)
  RkM <- matrix(0, 0, q) # Cache crossprod(Rk, M) for efficiency
  num_components <- 0

  for (kcomp in seq_len(K)) {
    lamk <- lambda.vec[kcomp]
    if (verbose) {
      message(sprintf("\nExtracting component %d/%d [penalty=%s, lambda=%.3g]",
                      kcomp, K, penalty, lamk))
    }

    # Precompute Cholesky for GPLS+Ridge case if needed
    chol_Q_lam_I <- NULL
    if (penalty == "ridge" && use.gpls) {
      Q_lam_I <- Q + lamk * Id
      # Use tryCatch for potential non-positive definite matrix during Cholesky
      chol_Q_lam_I <- tryCatch(chol(Q_lam_I), error = function(e) {
        warning(paste("Cholesky decomposition failed for Q + lambda*I at component", kcomp,
                      "; matrix might not be positive definite. Error:", e$message))
        NULL # Or handle differently, e.g., use solve() as fallback?
      })
      if (is.null(chol_Q_lam_I)) {
         # Fallback or stop? Let's stop for now if Cholesky fails.
         stop("Cholesky decomposition failed. Cannot proceed with GPLS + Ridge.")
      }
    }

    svd_M <- svd(M, nu = 1, nv = 1)
    v_k   <- svd_M$u  # dimension p
    u_k   <- svd_M$v  # dimension q


    # Check degenerate
    if (all(v_k == 0) || all(u_k == 0)) {
      if (verbose) message("Degenerate M => stopping.")
      break
    }

    diffU <- diffV <- Inf
    iter  <- 1
    while ((diffU > tol || diffV > tol) && iter <= maxiter) {
      old_u <- u_k
      old_v <- v_k

      # 1) update u_k
      if (use.gpls) {
        tmp <- crossprod(M, Q %*% v_k)  # M^T (Q v_k)
      } else {
        tmp <- crossprod(M, v_k)        # M^T v_k
      }
      denom_tmp <- sqrt(sum(tmp^2))
      if (denom_tmp > 1e-15) { # Use tolerance check
        u_k <- drop(tmp / denom_tmp)
      } else {
        u_k <- rep(0, q)
      }
      diffU <- sqrt(sum((u_k - old_u)^2))

      # 2) update v_k by penalized regression w/ partial "residual" r
      # NOTE: r calculation differs slightly between L1/L2 in some literature
      # Here, r is target for v_k: M %*% u_k (standard) or Q %*% (M %*% u_k) (GPLS)??
      # The objective implicitly involves Q norm if use.gpls=TRUE

      if (penalty == "l1") {
        # For L1, the subproblem is simpler: argmin_v 0.5||r-v||^2 + lamk ||v||_1
        # Where r = M %*% u_k (standard) or r = Q %*% M %*% u_k (GPLS)??
        # Let's stick to Allen's formulation: r = M %*% u_k and Q enters via norm.
        # So, r = M %*% u_k for both cases here.
        r <- M %*% u_k
        if (!nonneg) {
          v_k_new <- sign(r) * pmax(abs(r) - lamk, 0)
        } else {
          v_k_new <- pmax(r - lamk, 0)
        }
      } else if (penalty == "ridge") {
        if (nonneg) {
          warning("nonneg option is ignored for ridge penalties")
        }
         # For L2: objective is 0.5 * v^T Q v - v^T (Q M u) + 0.5 lamk v^T v
         # => derivative: Qv - Q M u + lamk v = 0 => (Q + lamk I) v = Q (M u)
         r <- M %*% u_k # Let r = M u

        if (!use.gpls) {
          # standard ridge, Q=I => (I + lamk I) v = I r => v = r / (1 + lamk)
          v_k_new <- r / (1 + lamk)
        } else {
          # GPLS + ridge => solve (Q + lamk I) v = Q r using precomputed Cholesky
          Qr <- Q %*% r
          # Solve Ax=b where A = Q + lamk I, b = Qr, using L Lt x = b
          # Step 1: Solve L y = b for y (forward solve)
          # Step 2: Solve Lt x = y for x (back solve)
          # Base R `solve` with Cholesky factor does this: solve(chol_A, b)
          # Alternatively: backsolve(R, forwardsolve(L, b)) where chol gives R=Lt
          v_k_new <- solve(chol_Q_lam_I, solve(t(chol_Q_lam_I), Qr))

          # Old direct solve:
          # v_k_new <- solve(Q + lamk * Id, Qr)
        }
      } else {
        # Should not happen due to match.arg, but defensive coding
        stop("Unrecognized penalty.")
      }

      # Now rescale v_k by 2-norm or Q-norm
      if (use.gpls) {
        # Q-norm: sqrt(v^T Q v)
        denom_v <- sqrt(drop(crossprod(v_k_new, Q %*% v_k_new)))
      } else {
        # 2-norm: sqrt(v^T v)
        denom_v <- sqrt(sum(v_k_new^2))
      }
      if (denom_v > 1e-15) {
        v_k <- v_k_new / denom_v
      } else {
        v_k <- rep(0, p)
      }
      diffV <- sqrt(sum((v_k - old_v)^2))

      if (verbose && iter %% 10 == 0) {
        message(sprintf("   iter=%d diffU=%.2e diffV=%.2e", iter, diffU, diffV))
      }
      iter <- iter + 1
    } # end while

    if (all(v_k == 0)) {
      if (verbose) message("v_k => all zero, stopping.")
      break
    }

    # Factor z_k
    if (use.gpls) {
      z_k <- X %*% (Q %*% v_k) # z = X Q v
    } else {
      z_k <- X %*% v_k        # z = X v
    }

    # store
    V[, kcomp] <- v_k
    U[, kcomp] <- u_k
    Z[, kcomp] <- z_k
    num_components <- num_components + 1

    # deflation (SIMPLS)
    denom_zk <- drop(crossprod(z_k, z_k))
    if (denom_zk < 1e-15) {
      if (verbose) message("degenerate z_k => stopping deflation.")
      break
    }

    r_k_new <- crossprod(X, z_k) / denom_zk
    Rk <- cbind(Rk, r_k_new)

    # Update RkM efficiently: calculate only the new column cross-product
    r_k_new_M <- crossprod(r_k_new, M)
    RkM <- rbind(RkM, r_k_new_M)

    # Calculate Rk^T Rk and its inverse using Cholesky
    RkTRk <- crossprod(Rk, Rk)
    chol_RkTRk <- tryCatch(chol(RkTRk), error = function(e) {
        warning(paste("Cholesky decomposition failed for Rk^T Rk at component", kcomp,
                      "; matrix might be ill-conditioned. Error:", e$message))
        NULL
    })

    if (is.null(chol_RkTRk)) {
       warning("Deflation unstable due to failure in Cholesky(Rk^T Rk), stopping.")
       break
    }

    # Use chol2inv for potentially better numerical stability / efficiency
    inv_RkTRk <- chol2inv(chol_RkTRk)

    # Old condition check (can still use rcond if needed, but chol gives direct check)
    # condRk <- rcond(RkTRk)
    # if (condRk < .Machine$double.eps) {
    #   warning("near-singular Rk^T Rk => deflation unstable, stopping.")
    #   break
    # }
    # inv_RkTRk <- solve(RkTRk) # Replaced by chol2inv

    # Deflate M using cached RkM
    M <- M - Rk %*% (inv_RkTRk %*% RkM) # M = M - Rk (Rk'Rk)^{-1} Rk' M

    # Recompute crossprod(Rk, M) after deflation so RkM reflects the
    # updated (deflated) M for the next iteration
    RkM <- crossprod(Rk, M)
  } # end for K

  if (num_components < K) {
    # trim
    V <- V[, seq_len(num_components), drop = FALSE]
    U <- U[, seq_len(num_components), drop = FALSE]
    Z <- Z[, seq_len(num_components), drop = FALSE]
  }

  list(V = V, U = U, Z = Z, num_components = num_components)
}

#' Regularised / Generalised Partial Least Squares (RPLS / GPLS)
#'
#' Implements the algorithm of Allen *et al.* (2013) for supervised
#' dimension-reduction with optional sparsity (\eqn{\ell_1}) or ridge
#' (\eqn{\ell_2}) penalties **and** the generalised extension that operates
#' in a user-supplied quadratic form \eqn{Q}.
#'
#' Unlike `genpls()`, which handles separate row and column metrics (`Mx`,
#' `Ax`, `My`, `Ay`) with a Gram–Schmidt orthogonalisation step, `rpls()`
#' uses a single metric `Q` and the simpler penalised updates of Allen et al.
#'
#' @section Method:
#' The routine follows Algorithm 1 of Allen *et al.* (2013, *Stat.
#' Anal. Data Min.*, 6 : 302–314) — see the paper for details. Briefly, with
#' \eqn{C = X^\top Y} the cross-product of the (preprocessed) blocks, each
#' component maximises
#' \deqn{\max_{u,v}\; v^\top Q C u - \lambda \, P(v)}
#' with \eqn{Q = I_p} for standard RPLS. The alternating updates are:
#' \eqn{u \leftarrow C^\top Q v / \|C^\top Q v\|_2}, then a penalised
#' (possibly non-negative) regression for \eqn{v}, normalised in the
#' \eqn{Q}-norm.
#'
#' @param X Numeric matrix \eqn{(n \times p)} — predictors.
#' @param Y Numeric matrix \eqn{(n \times q)} — responses.
#' @param K Integer, number of latent factors to extract. Default `2`.
#' @param lambda Scalar or length-\code{K} numeric vector of penalties.
#' @param penalty Either `"l1"` (lasso) or `"ridge"`.
#' @param Q Optional positive-(semi)definite \eqn{p \times p} matrix
#'   inducing *generalised* PLS. `NULL` means identity.
#' @param nonneg Logical, force non-negative loadings when
#'   \code{penalty = "l1"}. Note: This option is currently ignored when
#'   \code{penalty = "ridge"}.
#' @param tol Relative tolerance for the inner iterations convergence check.
#'   Default `1e-6`.
#' @param maxiter Maximum number of inner iterations per component. Default `200`.
#' @param verbose Logical; print progress messages during component extraction.
#'   Default `FALSE`.
#' @param preproc_x,preproc_y Optional \pkg{multivarious} preprocessing
#'   objects (see \code{\link[multivarious]{fit_transform}}). By default they pass
#'   the data through unchanged using \code{pass()}.
#' @param ... Further arguments (e.g., custom stopping criteria if implemented)
#'   are stored in the returned object (they are not used by
#'   \code{fit_rpls}).
#'
#' @return An object of class \code{c("rpls","cross_projector","projector")}
#'   with at least the elements
#'   \describe{
#'     \item{vx}{\eqn{p \times K} matrix of X-loadings.}
#'     \item{vy}{\eqn{q \times K} matrix of Y-loadings.}
#'     \item{ncomp}{Number of components actually extracted (may be < K).}
#'     \item{penalty}{Penalty type used (`"l1"` or `"ridge"`).}
#'     \item{lambda}{The `lambda` value(s) used.}
#'     \item{tol}{The convergence tolerance used.}
#'     \item{maxiter}{The maximum number of inner iterations used per component.}
#'     \item{nonneg}{Logical flag for the non-negativity constraint on `l1`.}
#'     \item{Q_used}{Logical; `TRUE` if a custom `Q` metric was supplied to
#'       induce generalised RPLS, `FALSE` for standard (identity-metric)
#'       RPLS. Note the field is named `Q_used`, not `Q`: the metric matrix
#'       itself is not retained on the returned object.}
#'     \item{verbose}{The `verbose` flag used.}
#'     \item{preproc_x, preproc_y}{Pre-processing transforms used.}
#'   }
#'
#'   `rpls` objects store no row (X-)scores: the factors \eqn{z_k} are
#'   computed internally during fitting to drive deflation and are then
#'   discarded, so use \code{multivarious::project()} on `X` to obtain
#'   scores for any rows of interest.
#'
#'   The object supports \code{predict()}, \code{project()},
#'   \code{transfer()}, \code{coef()} and other \pkg{multivarious} generics.
#'
#' @examples
#' # Generate sample data
#' set.seed(123)
#' n <- 50
#' p <- 20
#' q <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' Y <- X[, 1:5] %*% matrix(rnorm(5 * q), 5, q) + matrix(rnorm(n * q), n, q)
#'
#' # Fit regularized PLS with L1 penalty
#' fit_l1 <- rpls(X, Y, K = 3, lambda = 0.1, penalty = "l1")
#' print(fit_l1)
#'
#' # Fit regularized PLS with ridge penalty
#' fit_ridge <- rpls(X, Y, K = 3, lambda = 0.1, penalty = "ridge")
#' print(fit_ridge)
#'
#' @references Allen, G. I., Peterson, C., Vannucci, M., &
#'   Maletić-Savatić, M. (2013). *Regularized Partial Least Squares with
#'   an Application to NMR Spectroscopy.* **Statistical Analysis and Data
#'   Mining, 6(4)**, 302-314. DOI:10.1002/sam.11169.
#'
#' @importFrom multivarious cross_projector fit_transform pass
#' @export
rpls <- function(X, Y,
                 K         = 2,
                 lambda    = 0.1,
                 penalty   = c("l1", "ridge"),
                 Q         = NULL,
                 nonneg    = FALSE,
                 preproc_x = multivarious::pass(),
                 preproc_y = multivarious::pass(),
                 tol       = 1e-6,
                 maxiter   = 200,
                 verbose   = FALSE,
                 ...) {

  # Ensure penalty is valid early
  penalty  <- match.arg(penalty)

  # Preprocess X and Y using multivarious fit_transform API
  ft_x <- multivarious::fit_transform(preproc_x, as.matrix(X))
  ft_y <- multivarious::fit_transform(preproc_y, as.matrix(Y))
  proc_x <- ft_x$preproc
  proc_y <- ft_y$preproc
  Xp <- ft_x$transformed
  Yp <- ft_y$transformed

  # Call the internal fitting engine
  res_fit  <- fit_rpls(X = Xp, Y = Yp, K = K, lambda = lambda, penalty = penalty,
                       Q = Q, nonneg = nonneg,
                       tol = tol, maxiter = maxiter, verbose = verbose)

  # Package results into a cross_projector object
  out <- multivarious::cross_projector(
    vx        = res_fit$V,
    vy        = res_fit$U,
    preproc_x = proc_x,
    preproc_y = proc_y,
    classes   = "rpls",      # Add specific class first
    penalty   = penalty,
    ncomp     = res_fit$num_components,
    lambda    = lambda,      # Store lambda used
    Q_used    = !is.null(Q), # Indicate if Q was provided
    nonneg    = nonneg,      # Store nonneg flag
    tol       = tol,         # Store tolerance
    maxiter   = maxiter,     # Store max iterations
    verbose   = verbose,     # Store verbose flag
    ...                      # Store any other args passed
  )

  # Add back the generic classes if cross_projector doesn't handle inheritance fully
  # class(out) <- c("rpls", "cross_projector", "projector") # Ensure class hierarchy
  # Assuming cross_projector sets classes correctly (e.g., using c(classes, "cross_projector", "projector"))

  out
}
