#' Compatibility wrapper for multivarious::transfer
#'
#' Provides a transfer method that accepts `source`/`target` arguments
#' used in the tests. The method forwards to `multivarious`'s
#' transfer method which expects `from`/`to`.
#'
#' @param x A cross_projector object.
#' @param new_data New data to transfer.
#' @param from Source space ("X" or "Y").
#' @param to Target space ("X" or "Y").
#' @param opts Options list passed to multivarious::transfer.
#' @param ... Additional arguments. Legacy parameters `source` and `target`
#'   are accepted here for backwards compatibility but are deprecated;
#'   use `from` and `to` instead.
#' @return Matrix with transferred data.
#'
#' @examples
#' # Generate sample data
#' set.seed(123)
#' n <- 50
#' X <- matrix(rnorm(n * 10), n, 10)
#' Y <- matrix(rnorm(n * 8), n, 8)
#'
#' # Create a cross projector using rpls
#' fit <- rpls(X, Y, K = 2)
#'
#' # Transfer new X data to Y space
#' new_X <- matrix(rnorm(10 * 10), 10, 10)
#' transferred <- transfer(fit, new_X, from = "X", to = "Y")
#'
#' @importFrom multivarious transfer
#' @importFrom utils getS3method
#' @export
transfer.cross_projector <- function(x, new_data, from = NULL, to = NULL,
                                     opts = list(), ...) {

  # Handle legacy parameter names (source/target -> from/to) from ...
  dots <- list(...)
  if (!is.null(dots$source) && is.null(from)) {
    from <- dots$source
  }
  if (!is.null(dots$target) && is.null(to)) {
    to <- dots$target
  }

  # Set defaults if not provided
  if (is.null(from)) from <- "X"
  if (is.null(to)) to <- "Y"

  # Match arguments
  from <- match.arg(from, c("X", "Y"))
  to   <- match.arg(to, c("X", "Y"))

  # Call multivarious transfer method directly to avoid method dispatch recursion
  getS3method("transfer", "cross_projector", envir = asNamespace("multivarious"))(
    x, new_data, from = from, to = to, opts = opts
  )
}
