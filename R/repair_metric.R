#' Repair a metric matrix explicitly
#'
#' Returns a positive (semi)definite version of `A` together with a
#' diagnostic report of what was done. This is the explicit counterpart of
#' the `constraints_remedy` argument of [genpca()]: nothing in the package
#' repairs a metric silently, and this function lets you inspect the repair
#' before using the result.
#'
#' @param A A square symmetric matrix (base matrix or `Matrix`). Asymmetry
#'   beyond roundoff is an error (see the relative asymmetry test in
#'   `symmetrize_or_stop()`); it is not something a PSD repair should hide.
#' @param method `"ridge"` adds a diagonal loading that makes the
#'   matrix positive definite (Gershgorin-based shift, falling back to
#'   `Matrix::nearPD()` for small dense matrices); `"clip"` projects onto the
#'   PSD cone by zeroing negative eigenvalues (dense eigendecomposition; refuses
#'   large sparse input); `"identity"` replaces an indefinite matrix by the
#'   identity (the report still describes the input).
#' @param rtol Relative tolerance: eigenvalues above `-rtol * scale(A)` count
#'   as non-negative for `"ridge"` and `"identity"`, which then return the
#'   matrix unchanged. `"clip"` always removes negative eigenvalues, regardless
#'   of this tolerance (up to reconstruction roundoff). Default
#'   `sqrt(.Machine$double.eps)`.
#' @param name Label used in messages.
#' @param diag_maxn Largest dimension for which the report computes the full
#'   spectrum (minimum eigenvalue, rank, condition number); above it an
#'   iterative minimum eigenvalue estimate and a Gershgorin bound are reported,
#'   with rank and condition number unavailable.
#' @return The repaired matrix (a `Matrix`), with attribute `"repair_report"`
#'   of class `"metric_repair_report"`: a list with `name`, `method`,
#'   `changed`, `n`, `min_eigenvalue_before`, `min_eigenvalue_after`,
#'   `gershgorin_bound_before`, `shift` (diagonal loading added by `"ridge"`,
#'   `NA` for `"clip"`), `rank`, `condition_number` and `rtol`.
#' @examples
#' A <- matrix(c(1, 2, 2, 1), 2)           # eigenvalues 3 and -1
#' B <- repair_metric(A, method = "ridge")
#' attr(B, "repair_report")
#' C <- repair_metric(A, method = "clip")
#' eigen(as.matrix(C))$values
#' @seealso [genpca()] (argument `constraints_remedy`)
#' @export
repair_metric <- function(A, method = c("ridge", "clip", "identity"),
                          rtol = .metric_rtol_default(), name = "A",
                          diag_maxn = 2000L) {
  method <- match.arg(method)
  A <- symmetrize_or_stop(A, name = name)
  n <- nrow(A)
  vals <- if (methods::is(A, "sparseMatrix")) A@x else as.numeric(A)
  if (any(!is.finite(vals))) stop("Metric contains non-finite values.", call. = FALSE)

  spectrum <- function(M) {
    if (Matrix::isDiagonal(M)) {
      ev <- as.numeric(Matrix::diag(M))
      pos <- ev[ev > rtol * max(ev, 0)]
      return(list(min = min(ev), rank = length(pos),
                  cond = if (length(pos) && min(pos) > 0) max(pos) / min(pos) else NA_real_))
    }
    if (n <= diag_maxn) {
      ev <- eigen(as.matrix(M), symmetric = TRUE, only.values = TRUE)$values
      pos <- ev[ev > rtol * max(ev, 0)]
      list(min = min(ev), rank = length(pos),
           cond = if (length(pos) && min(pos) > 0) max(pos) / min(pos) else NA_real_)
    } else {
      mn <- tryCatch(.top_eigs_sym(M, 1, "SA")$values[1], error = function(e) NA_real_)
      list(min = mn, rank = NA_integer_, cond = NA_real_)
    }
  }
  d <- as.numeric(Matrix::diag(A))
  gersh <- suppressWarnings(min(d - (Matrix::rowSums(abs(A)) - abs(d))))
  before <- spectrum(A)
  changed <- !is_psd(A, rtol = rtol)

  if (method == "clip") {
    B <- clip_psd(A, name = name)
    changed <- any(A != B)
  } else if (changed) {
    B <- switch(method,
                clip = clip_psd(A, name = name),
                identity = Matrix::Diagonal(n),
                ridge = ensure_spd(A, tol = rtol, name = name))
  } else {
    B <- A
  }
  after <- spectrum(B)
  shift <- if (!changed) 0 else if (method == "ridge") {
    as.numeric(Matrix::diag(B)[1] - d[1])
  } else NA_real_

  report <- structure(list(
    name = name, method = method, changed = changed, n = n,
    min_eigenvalue_before = before$min, min_eigenvalue_after = after$min,
    gershgorin_bound_before = gersh, shift = shift,
    rank = after$rank, condition_number = after$cond, rtol = rtol
  ), class = "metric_repair_report")
  attr(B, "repair_report") <- report
  B
}

#' @export
print.metric_repair_report <- function(x, ...) {
  cat("Metric repair report for ", x$name, " (", x$n, " x ", x$n, ")\n", sep = "")
  cat("  method:                ", x$method, "\n", sep = "")
  cat("  changed:               ", x$changed, "\n", sep = "")
  cat("  min eigenvalue before: ", format(x$min_eigenvalue_before, digits = 4), "\n", sep = "")
  cat("  min eigenvalue after:  ", format(x$min_eigenvalue_after, digits = 4), "\n", sep = "")
  cat("  Gershgorin bound:      ", format(x$gershgorin_bound_before, digits = 4), "\n", sep = "")
  if (x$method == "ridge") cat("  diagonal shift:        ", format(x$shift, digits = 4), "\n", sep = "")
  cat("  rank:                  ", x$rank, "\n", sep = "")
  cat("  condition number:      ", format(x$condition_number, digits = 4), "\n", sep = "")
  cat("  relative tolerance:    ", format(x$rtol, digits = 3), "\n", sep = "")
  invisible(x)
}

#' Signal that a metric was repaired
#'
#' Emits a warning of class `genpca_metric_repaired` carrying the report, so
#' callers can catch it (`withCallingHandlers`) or silence it selectively.
#' @keywords internal
#' @noRd
.warn_metric_repaired <- function(report) {
  msg <- if (report$method == "identity") {
    sprintf("Matrix %s is not positive semi-definite (min eigenvalue %s); replaced by the identity as requested by constraints_remedy = \"identity\".",
            report$name, format(report$min_eigenvalue_before, digits = 3))
  } else {
    sprintf("Matrix %s is not positive semi-definite (min eigenvalue %s); repaired with constraints_remedy = \"%s\"%s. Inspect with repair_metric().",
            report$name, format(report$min_eigenvalue_before, digits = 3), report$method,
            if (report$method == "ridge" && is.finite(report$shift)) sprintf(" (diagonal shift %s)", format(report$shift, digits = 3)) else "")
  }
  warning(structure(class = c("genpca_metric_repaired", "warning", "condition"),
                    list(message = msg, call = NULL, report = report)))
}
