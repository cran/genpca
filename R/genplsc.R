#' Canonical Generalized PLS (alias)
#'
#' Convenience alias for `genpls()`; computes canonical generalized PLS
#' (PLS-SVD/GPLSSVD). See `?genpls` for full documentation.
#'
#' @inheritParams genpls
#' @return An object of class `c("genpls", "cross_projector", "projector")`
#'   with the same structure as `genpls()` returns (X-/Y-weights `vx`/`vy`,
#'   singular values `d`, generalized weights `p`/`q`, scores `fi`/`fj`,
#'   latent variables `lx`/`ly`, `ncomp`, and `backend`); see `?genpls` for
#'   the definition of each slot.
#' @seealso [genpls()]
#' @references
#' Beaton, D. (2020). Generalized eigen, singular value, and partial least
#' squares decompositions: The GSVD package. (Eqs. 10-14). arXiv:2010.14734.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(60 * 5), 60, 5)
#' Y <- matrix(rnorm(60 * 4), 60, 4)
#' fit <- genplsc(X, Y, ncomp = 2,
#'                preproc_x = multivarious::center(),
#'                preproc_y = multivarious::center())
#' fit$d
#' @export
genplsc <- function(X, Y,
                    Ax = NULL, Ay = NULL,
                    Mx = NULL, My = NULL,
                    ncomp = 2,
                    preproc_x = multivarious::pass(),
                    preproc_y = multivarious::pass(),
                    svd_backend = c("eigencore", "irlba", "RSpectra"),
                    svd_opts = list(tol = 1e-7, maxitr = 1000),
                    constraints_remedy = c("error", "ridge", "clip", "identity"),
                    verbose = FALSE) {
  genpls(X = X, Y = Y,
         Ax = Ax, Ay = Ay,
         Mx = Mx, My = My,
         ncomp = ncomp,
         preproc_x = preproc_x,
         preproc_y = preproc_y,
         svd_backend = svd_backend,
         svd_opts = svd_opts,
         constraints_remedy = constraints_remedy,
         verbose = verbose)
}
