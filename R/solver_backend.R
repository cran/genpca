#' Iterative eigen/SVD backend
#'
#' Every partial eigen and singular value computation in genpca goes through
#' `.top_eigs_sym()` and `.top_svd()`, which wrap the **eigencore** package.
#' Keeping the solver in one place fixes four behaviours that differ between
#' iterative solvers: tolerance/iteration options are passed explicitly,
#' function operators are wrapped as `eigencore::linear_operator()` with
#' block (GEMM-style) callbacks, the caller's `.Random.seed` is never
#' disturbed, and incomplete results raise a `genpca_solver_nonconvergence`
#' error based on the returned `nconv` (the
#' certificate's `passed` flag is withheld for matrix-free operators, whose
#' norm bound is only estimated).
#'
#' @keywords internal
#' @noRd
NULL

.with_rng_guard <- function(expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  expr
}

# Matrices go to the solver as base double matrices or dgCMatrix.
.solver_matrix <- function(A) {
  if (methods::is(A, "sparseMatrix")) return(as_dgc(A))
  if (inherits(A, "Matrix")) return(as.matrix(A))
  if (!is.double(A)) storage.mode(A) <- "double"
  A
}

# eigencore calls operator callbacks as apply(X, alpha, beta, Y) with X an
# n x b block and expects a double matrix back. genpca's operator closures
# take (x, args) and may return numeric vectors, so reshape here.
# force() matters: callers rebind the variable holding `f` after this
# closure is created.
.gemm_callback <- function(f, nrow_out, args = NULL) {
  force(f); force(nrow_out); force(args)
  function(X, alpha = 1, beta = 0, Y = NULL, ...) {
    X <- as.matrix(X)
    out <- matrix(as.numeric(f(X, args)), nrow = nrow_out, ncol = ncol(X))
    if (!identical(alpha, 1)) out <- alpha * out
    if (!is.null(Y) && !identical(beta, 0)) out <- out + beta * as.matrix(Y)
    out
  }
}

#' Top-k eigenpairs of a symmetric matrix
#' @param A symmetric matrix (base, dense Matrix or sparse Matrix)
#' @param k number of eigenpairs
#' @param which "LA" (largest algebraic), "LM" (largest magnitude) or "SA"
#'   (smallest algebraic)
#' @param tol convergence tolerance
#' @param maxit optional iteration limit
#' @return list(values, vectors, nconv, converged)
#' @noRd
.top_eigs_sym <- function(A, k, which = c("LA", "LM", "SA"), tol = 1e-8, maxit = NULL) {
  which <- match.arg(which)
  target <- switch(which,
                   LA = eigencore::largest(),
                   LM = eigencore::largest_magnitude(),
                   SA = eigencore::smallest())
  A <- .solver_matrix(A)
  k <- as.integer(k)
  P <- eigencore::eigen_problem(A, structure = eigencore::hermitian(), target = target)
  fit <- .with_rng_guard(solve(P, k = k, tol = tol, maxit = maxit))
  .require_solver_convergence(fit$nconv, min(k, nrow(A)), "eigen")
  list(values = as.numeric(fit$values),
       vectors = fit$vectors,
       nconv = fit$nconv,
       converged = isTRUE(fit$nconv >= min(k, nrow(A))))
}

#' Top-k singular triplets of a matrix or operator
#' @param A matrix, or a function `(x, args)` computing `A %*% x`
#' @param k number of triplets
#' @param nu,nv number of left/right vectors (0 to skip)
#' @param tol convergence tolerance
#' @param adjoint for function `A`: function `(x, args)` computing `t(A) %*% x`
#' @param dim for function `A`: `c(nrow, ncol)`
#' @param args extra argument passed to the callbacks
#' @return list(d, u, v, nconv, converged)
#' @noRd
.top_svd <- function(A, k, nu = k, nv = k, tol = 1e-8, adjoint = NULL,
                     dim = NULL, args = NULL) {
  if (is.function(A)) {
    stopifnot(is.function(adjoint), length(dim) == 2L)
    dim <- as.integer(dim)
    A <- eigencore::linear_operator(
      dim = dim,
      apply = .gemm_callback(A, dim[1], args),
      apply_adjoint = .gemm_callback(adjoint, dim[2], args))
  } else {
    A <- .solver_matrix(A)
  }
  vectors <- if (nu > 0 && nv > 0) "both" else if (nu > 0) "left" else if (nv > 0) "right" else "none"
  k <- as.integer(k)
  fit <- .with_rng_guard(eigencore::svd_partial(A, rank = k, target = eigencore::largest(),
                                                tol = tol, vectors = vectors))
  .require_solver_convergence(fit$nconv, k, "SVD")
  list(d = as.numeric(fit$d), u = fit$u, v = fit$v, nconv = fit$nconv,
       converged = isTRUE(fit$nconv >= k))
}

# Refuse incomplete results centrally. Callers with a dense fallback already
# catch this error; matrix-free callers propagate it without allocating a
# potentially unbounded dense matrix.
.require_solver_convergence <- function(nconv, k, solver) {
  if (length(nconv) != 1L || !is.finite(nconv) || nconv < k) {
    stop(structure(list(
      message = paste0("eigencore ", solver, " did not converge: ",
                       if (length(nconv) == 1L) nconv else "unknown",
                       " of ", k, " requested components converged. ",
                       "Try a less stringent solver tolerance."),
      call = NULL, nconv = nconv, requested = k),
      class = c("genpca_solver_nonconvergence", "error", "condition")))
  }
  invisible(TRUE)
}
