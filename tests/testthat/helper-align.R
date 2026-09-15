# Align signs of A_test's columns to A_ref's (flip columns whose dot product
# with the reference column is negative). Shared by several test files.
align_signs <- function(A_ref, A_test) {
  if (is.null(A_ref) || is.null(A_test)) return(A_test)
  A_ref <- as.matrix(A_ref)
  A_test <- as.matrix(A_test)
  if (ncol(A_ref) == 0 || ncol(A_test) == 0) return(A_test)
  stopifnot(ncol(A_ref) == ncol(A_test))
  s <- sign(colSums(A_ref * A_test))
  s[s == 0] <- 1
  A_test %*% diag(s, ncol(A_ref))
}

# Align columns of A to B up to permutation and sign; return correlation diag.
align_perm_sign <- function(A, B) {
  stopifnot(ncol(A) == ncol(B))

  # Check if clue package is available
  if (!requireNamespace("clue", quietly = TRUE)) {
    warning("Package 'clue' not available for optimal alignment. Using simpler matching.")
    # Simple greedy matching based on absolute correlations
    C <- abs(cor(A, B))
    perm <- integer(ncol(A))
    for (i in seq_len(ncol(A))) {
      j <- which.max(C[i, ])
      perm[i] <- j
      C[, j] <- 0  # Mark as used
    }
    return(list(corr = diag(abs(cor(A, B[, perm]))), perm = perm))
  }

  # Optimal matching using Hungarian algorithm
  C <- abs(cor(A, B))
  perm <- clue::solve_LSAP(C, maximum = TRUE)
  list(corr = C[cbind(seq_len(ncol(A)), perm)],
       perm = perm)
}
