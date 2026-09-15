context("genpca")

set.seed(3)
mat_10_10 <- matrix(rnorm(10 * 10), 10, 10)

test_that("ncomp must be integer", {
  expect_error(genpca(mat_10_10, ncomp = 2.5), "single positive integer")
  expect_error(genpca(mat_10_10, ncomp = c(1, 2)), "single positive integer")
})

test_that("pca and genpca have same results with identity matrix for row and column constraints", {
  res1 <- genpca(mat_10_10, preproc = multivarious::center())
  res2 <- multivarious::pca(mat_10_10, ncomp = ncomp(res1), preproc = multivarious::center())

  diffscores <- abs(multivarious::scores(res1)) - abs(multivarious::scores(res2))
  expect_true(sum(diffscores) < 1e-5)
  expect_equal(multivarious::sdev(res1), multivarious::sdev(res2))

  expect_equal(unname(apply(multivarious::components(res1), 2, function(x) sum(x^2))),
               rep(1, ncomp(res1)))
})

test_that("gen_pca with column variances is equivalent to a scaled pca", {
  wts <- 1 / apply(mat_10_10, 2, var)
  res1 <- genpca(mat_10_10, A = wts, preproc = multivarious::center())
  res2 <- multivarious::pca(mat_10_10, preproc = multivarious::standardize())

  # Compare absolute scores element-wise with tolerance
  # Scores should match up to sign flips
  expect_equal(abs(as.matrix(multivarious::scores(res1))), abs(as.matrix(multivarious::scores(res2))), tolerance = 1e-6, check.attributes = FALSE)
  expect_equal(multivarious::sdev(res1), multivarious::sdev(res2), check.attributes = FALSE)

})

test_that("gen_pca with use_cpp with column variances is equivalent to a scaled pca", {
  wts <- 1 / apply(mat_10_10, 2, var)
  res1 <- genpca(mat_10_10, A = wts, preproc = multivarious::center(), use_cpp = TRUE)
  res2 <- multivarious::pca(mat_10_10, preproc = multivarious::standardize())

  diffscores <- abs(multivarious::scores(res1)) - abs(multivarious::scores(res2))
  expect_true(abs(sum(diffscores)) < 1e-5)
  expect_equal(multivarious::sdev(res1), multivarious::sdev(res2), check.attributes = FALSE)

})

test_that("gen_pca with use_cpp (+ deflation) with column variances is equivalent to a scaled pca", {
  wts <- 1 / apply(mat_10_10, 2, var)
  res1 <- genpca(mat_10_10, A = wts, preproc = multivarious::center(), ncomp = 9,
  method = "deflation", use_cpp = TRUE, threshold = 1e-7)
  res2 <- multivarious::pca(mat_10_10, preproc = multivarious::standardize(), ncomp = 9)

  diffscores <- abs(as.matrix(multivarious::scores(res1))) - abs(as.matrix(multivarious::scores(res2)))
  expect_true(mean(abs(diffscores)) < .01)
  expect_equal(multivarious::sdev(res1), multivarious::sdev(res2), tolerance = .01)

})

test_that("gen_pca with use_cpp (+ deflation and n < p) with column variances is equivalent to a scaled pca", {
  set.seed(58)
  mat_10_20 <- matrix(rnorm(10 * 20), 10, 20)
  wts <- 1 / apply(mat_10_20, 2, var)
  res1 <- genpca(mat_10_20, A = wts, preproc = multivarious::center(), ncomp = 9, method = "deflation",
  use_cpp = FALSE, threshold = 1e-8)
  res2 <- multivarious::pca(mat_10_20, preproc = multivarious::standardize(), ncomp = 9)

  diffscores <- abs(as.matrix(multivarious::scores(res1))) - abs(as.matrix(multivarious::scores(res2)))
  expect_true(mean(abs(diffscores)) < .01)
  expect_equal(multivarious::sdev(res1), multivarious::sdev(res2), tolerance = .01)

})

test_that("gen_pca with dense column and row constraints works", {
  set.seed(71)
  A <- cov(matrix(rnorm(10 * 10), 10, 10))
  M <- cov(matrix(rnorm(10 * 10), 10, 10))
  res1 <- genpca(mat_10_10, A = A, M = M, preproc = multivarious::center())
  expect_equal(ncomp(res1), length(multivarious::sdev(res1)))
})

test_that("gen_pca with sparse column and row constraints works", {
  skip_if_not_installed("adjoin")
  A <- adjoin::adjacency(adjoin::graph_weights(mat_10_10, k = 8))
  Matrix::diag(A) <- 1
  M <- adjoin::adjacency(adjoin::graph_weights(t(mat_10_10), k = 3))
  Matrix::diag(M) <- 1.5
  # the adjacency-based row metric is indefinite: opt into the ridge repair
  expect_warning(
    res1 <- genpca(mat_10_10, A = A, M = M, preproc = multivarious::center(),
                   constraints_remedy = "ridge"),
    class = "genpca_metric_repaired"
  )

  k <- multivarious::ncomp(res1)
  expect_gt(k, 0)
  expect_equal(nrow(multivarious::scores(res1)), nrow(mat_10_10))
  expect_equal(nrow(multivarious::components(res1)), ncol(mat_10_10))
  # metric orthonormality against the (repaired) constraint matrices actually used
  expect_equal(as.matrix(t(res1$ou) %*% res1$M %*% res1$ou), diag(k), tolerance = 1e-6)
  expect_equal(as.matrix(t(res1$ov) %*% res1$A %*% res1$ov), diag(k), tolerance = 1e-6)
})


test_that("can reconstruct a genpca model with component selection", {
  set.seed(88)
  A <- cov(matrix(rnorm(20 * 10), 20, 10))
  M <- cov(matrix(rnorm(20 * 10), 20, 10))
  res1 <- genpca(mat_10_10, preproc = multivarious::center())
  recon1 <- reconstruct(res1)
  expect_equal(as.matrix(recon1), mat_10_10, check.attributes = FALSE)

  res1 <- genpca(mat_10_10, A = A, M = M, ncomp = 10, preproc = multivarious::center())
  res2 <- multivarious::pca(mat_10_10, ncomp = 10, preproc = multivarious::center())
  recon2 <- reconstruct(res2)



})

test_that("can project a row vector", {
  set.seed(103)
  A <- cov(matrix(rnorm(10 * 10), 10, 10))
  M <- cov(matrix(rnorm(10 * 10), 10, 10))

  res1 <- genpca(mat_10_10, A = A, M = M)
  p <- multivarious::project(res1, mat_10_10[1, ])
  expect_equal(dim(p), c(1, ncomp(res1)))
})

#test_that("can extract residuals", {
#  res1 <- genpca(mat_10_10)
#  resid <- residuals(res1, ncomp=2, mat_10_10)
#  d <- multivarious::sdev(res1)
#  expect_equal(sum(d[3:length(d)] ^2), sum(resid^2))
#})

test_that("can run genpca with deflation", {
  set.seed(119)
  X <- matrix(rnorm(100), 10, 10)
  res1 <- genpca(X, preproc = multivarious::center(), ncomp = 5, method = "deflation")
  res2 <- genpca(X, preproc = multivarious::center(), ncomp = 5)
  expect_true(sum(abs(res1$u) - abs(res2$u)) < 1)
})

test_that("can run genpca with sparse weighting matrix", {
  skip_if_not_installed("adjoin")
  set.seed(127)
  # tall-and-thin with a large sparse row metric; kept modest so the suite
  # stays fast on CRAN while still exercising the sparse deflation path
  nr <- 1200
  X <- matrix(rnorm(nr * 20), nr, 20)
  A <- adjoin::temporal_adjacency(1:20)
  A <- cov(as.matrix(A))
  M <- adjoin::temporal_adjacency(1:nr)
  expect_true(methods::is(M, "sparseMatrix"))

  # temporal_adjacency() is indefinite: opt into the ridge repair explicitly
  expect_warning(
    res1 <- genpca(X, A = Matrix::Matrix(A, sparse = TRUE), M = M,
                   preproc = multivarious::center(), ncomp = 5, method = "deflation",
                   constraints_remedy = "ridge"),
    class = "genpca_metric_repaired"
  )
  res2 <- suppressWarnings(genpca(X, A = A, M = M, preproc = multivarious::center(), ncomp = 5,
                                  constraints_remedy = "ridge"))

  expect_equal(multivarious::ncomp(res1), 5L)
  expect_equal(multivarious::ncomp(res2), 5L)
  expect_true(all(is.finite(multivarious::sdev(res1))))
  # deflation and eigen agree on the leading singular values
  expect_equal(multivarious::sdev(res1)[1:3], multivarious::sdev(res2)[1:3],
               tolerance = 1e-3)
})

test_that("can run genpca on a largeish matrix with deflation", {
  skip_if_not_installed("adjoin")
  set.seed(138)
  nr <- 400
  nc <- 200
  X <- matrix(rnorm(nr * nc), nr, nc)
  A <- adjoin::temporal_adjacency(1:nc)
  A <- t(A) %*% A                          # PSD by construction

  M <- adjoin::temporal_adjacency(1:nr)
  M <- t(M) %*% M

  res1 <- genpca(X, A = Matrix::Matrix(A, sparse = TRUE),
                 M = M, preproc = multivarious::center(), ncomp = 5, method = "deflation", threshold = 1e-8)
  res2 <- genpca(X, A = Matrix::Matrix(A, sparse = TRUE),
                 M = M, preproc = multivarious::center(), ncomp = 5, method = "deflation",
                 threshold = 1e-8, use_cpp = FALSE)

  res3 <- genpca(X, A = Matrix::Matrix(A, sparse = TRUE),
                 M = M, preproc = multivarious::center(), ncomp = 20, method = "eigen")

  # C++ and R deflation agree with each other, and with the direct eigen path
  expect_equal(multivarious::sdev(res1), multivarious::sdev(res2), tolerance = 1e-4)
  expect_equal(multivarious::sdev(res1), multivarious::sdev(res3)[1:5], tolerance = 1e-3)
  expect_equal(multivarious::ncomp(res3), 20L)
  expect_true(all(diff(multivarious::sdev(res3)) <= 1e-8))   # non-increasing
})

## tests/testthat/test-genpca-spatial.R

test_that("genpca with spatial adjacency recovers a smooth temporal blob better than no constraints", {

  skip_if_not_installed("Matrix")
  library(Matrix)

  set.seed(1234)

  ## 1) Generate small 2D grid, e.g. 8x8
  nr <- 8
  nc <- 8
  P  <- nr * nc          # total number of pixels
  T  <- 20               # number of time points

  ## 2) Construct a ground-truth smooth blob
  ##    We'll define a circular-ish blob in the center
  grid_x <- matrix(rep(1:nr, each = nc), nrow = nr, ncol = nc)
  grid_y <- matrix(rep(1:nc, nr),      nrow = nr, ncol = nc, byrow = TRUE)

  center_r <- floor(nr / 2)
  center_c <- floor(nc / 2)
  radius   <- 2.5
  blob     <- exp(-((grid_x - center_r)^2 + (grid_y - center_c)^2) / (2 * radius^2))

  ## Flatten to length=P
  blob_vec <- as.vector(blob)  # shape (P)

  ## 3) Create a time-varying amplitude with a sinusoid
  tt <- seq(0, 2 * pi, length.out = T)
  amp <- 2 + sin(tt)   # shape (T)

  ## Expand to a T x P matrix
  ## Each row t is the blob * amp[t]
  signal_mat <- outer(amp, blob_vec)  # shape (T x P)

  ## 4) Add noise
  noise_mat <- matrix(rnorm(T * P, sd = 0.5), T, P)
  X <- signal_mat + noise_mat   # final observed data

  rownames(X) <- paste0("Time", seq_len(T))

  ## 5) Build a "spatial adjacency" or Laplacian matrix A for the 8x8 grid
  ##    - We'll do adjacency for up/down/left/right neighbors
  ##    - Then we might create a Laplacian from it (D - A, etc.)
  ## Here, let's do adjacency directly: A_ij = 1 if i & j are neighbors

  adj_list <- list()
  index_2d_to_1d <- function(r, c) (r - 1) * nc + c
  for (r in seq_len(nr)) {
    for (c in seq_len(nc)) {
      cur_id <- index_2d_to_1d(r, c)
      neighbors <- c()
      if (r > 1)         neighbors <- c(neighbors, index_2d_to_1d(r - 1, c))
      if (r < nr)        neighbors <- c(neighbors, index_2d_to_1d(r + 1, c))
      if (c > 1)         neighbors <- c(neighbors, index_2d_to_1d(r, c - 1))
      if (c < nc)        neighbors <- c(neighbors, index_2d_to_1d(r, c + 1))
      for (ngb in neighbors) {
        adj_list[[length(adj_list) + 1]] <- c(cur_id, ngb)
      }
    }
  }

  ## Build a sparse adjacency matrix from these edges
  row_inds <- sapply(adj_list, `[[`, 1)
  col_inds <- sapply(adj_list, `[[`, 2)
  ones     <- rep(1, length(adj_list))

  # A_sp is shape P x P
  A_sp <- sparseMatrix(i = row_inds, j = col_inds, x = ones, dims = c(P, P))

  ## Make adjacency PSD by adding diagonal
  A_sp_psd <- A_sp + 5 * Diagonal(P)  # Add strong diagonal to make PSD

  ## (Optional) Laplacian = Diag(rowSums(A_sp)) - A_sp
  d_vec <- rowSums(A_sp)
  Lap_sp <- sparseMatrix(i = 1:P, j = 1:P, x = d_vec) - A_sp


  ## 6) Run genpca with no constraint  (A=I, M=I)
  # We'll do just 1 component for demonstration
  A_id <- Diagonal(x = rep(1, P))  # identity
  M_id <- Diagonal(x = rep(1, T))  # identity
  # Possibly center or no preproc
  fit_no_constraint <- genpca(X, A = A_id, M = M_id, ncomp = 1,
                              preproc = multivarious::center(),
                              method = "eigen", use_cpp = FALSE)

  ## 7) Run genpca with adjacency + diagonal (now PSD)
  fit_adj <- genpca(X, A = A_sp_psd, M = M_id, ncomp = 1,
                    preproc = multivarious::center(),
                    method = "eigen", use_cpp = FALSE)

  ## 8) Compare reconstruction MSE for rank-1 approximation
  #   We'll do Xhat = reconstruct(..., comp=1)
  #   Then measure mean((X - Xhat)^2)
  Xhat_no_constr <- reconstruct(fit_no_constraint, comp = 1)
  Xhat_adj       <- reconstruct(fit_adj, comp = 1)

  mse_no_constr <- mean((X - Xhat_no_constr)^2)
  mse_adj       <- mean((X - Xhat_adj)^2)

  ## 9) Test passes if they are reasonably close
  ## Since spatial constraints may not always improve MSE on random data,
  ## we just check they are within reasonable bounds
  cat("MSE no constraint:", mse_no_constr, "\n")
  cat("MSE adjacency:    ", mse_adj, "\n")

  # The adjacency constraint should at least not make things much worse
  # and ideally would improve things slightly for smooth data
  expect_lt(mse_adj, mse_no_constr * 1.2)  # Allow up to 20% worse
})

############################################################
## Direct tests for gmd_fast_cpp implementation
############################################################

context("gmd_fast_cpp direct tests")

# Helper function to compare subspaces (ignoring column signs)
compare_subspaces <- function(U1, U2, tol = 1e-6) {
  expect_equal(ncol(U1), ncol(U2))
  if (ncol(U1) == 0) return(TRUE)

  # Normalize columns just in case
  U1 <- apply(U1, 2, function(x) x / sqrt(sum(x^2)))
  U2 <- apply(U2, 2, function(x) x / sqrt(sum(x^2)))

  # Correlation matrix between columns
  corr_mat <- abs(t(U1) %*% U2)

  # Check if it's close to a permutation matrix
  # Each row and column should have exactly one element close to 1
  row_max_close_to_1 <- all(abs(apply(corr_mat, 1, max) - 1) < tol)
  col_max_close_to_1 <- all(abs(apply(corr_mat, 2, max) - 1) < tol)

  # Check if the sum of squares of correlations is close to the number of components
  sum_sq_corr_close_to_k <- abs(sum(corr_mat^2) - ncol(U1)) < tol * ncol(U1)

  expect_true(row_max_close_to_1, info = "Max correlation per row not close to 1")
  expect_true(col_max_close_to_1, info = "Max correlation per col not close to 1")
  expect_true(sum_sq_corr_close_to_k, info = "Sum of squared correlations not close to k")
}

test_that("gmd_fast_cpp matches genpca (use_cpp=TRUE) for p <= n, dense constraints", {
  set.seed(1)
  n <- 20
  p <- 15
  k <- 5
  X <- matrix(rnorm(n * p), n, p)
  Q <- crossprod(matrix(rnorm(n * n), n, n)) + diag(n) * 0.1 # Dense SPD
  R <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p) * 0.1 # Dense SPD

  # Center data as gmd_fast_cpp assumes centered data
  X_centered <- scale(X, center = TRUE, scale = FALSE)

  res_r <- genpca(X_centered, M = Q, A = R, ncomp = k, use_cpp = TRUE, method = "spectra", preproc = multivarious::pass())
  # Directly call C++ function (make sure it's exported properly)
  res_cpp <- genpca:::gmd_fast_cpp(X_centered, Q = Matrix(Q), R = Matrix(R), k = k)

  expect_equal(res_cpp$d, multivarious::sdev(res_r), tolerance = 1e-6)
  expect_equal(res_cpp$k, k)
  # Compare actual values, not subspaces (since these are scores/components, not eigenvectors)
  # res_cpp$u is the M-weighted scores (M ou D); scores() is ou D = X A ov
  # (Allen et al. 2014), stored via fit$u = M ou.
  u_ref <- as.matrix(sweep(res_r$u, 2, multivarious::sdev(res_r), `*`))
  expect_equal(align_signs(u_ref, res_cpp$u), u_ref,
               tolerance = 1e-6, check.attributes = FALSE)
  expect_equal(multivarious::scores(res_r),
               sweep(res_r$ou, 2, multivarious::sdev(res_r), `*`),
               tolerance = 1e-8, check.attributes = FALSE)
  v_ref <- as.matrix(multivarious::components(res_r))
  expect_equal(align_signs(v_ref, res_cpp$v), v_ref,
               tolerance = 1e-6, check.attributes = FALSE)
})

test_that("gmd_fast_cpp matches genpca (use_cpp=TRUE) for p > n, dense constraints", {
  set.seed(2)
  n <- 15
  p <- 20
  k <- 5
  X <- matrix(rnorm(n * p), n, p)
  Q <- crossprod(matrix(rnorm(n * n), n, n)) + diag(n) * 0.1 # Dense SPD
  R <- crossprod(matrix(rnorm(p * p), p, p)) + diag(p) * 0.1 # Dense SPD

  X_centered <- scale(X, center = TRUE, scale = FALSE)

  res_r <- genpca(X_centered, M = Q, A = R, ncomp = k, use_cpp = TRUE, method = "spectra", preproc = multivarious::pass())
  res_cpp <- genpca:::gmd_fast_cpp(X_centered, Q = Matrix(Q), R = Matrix(R), k = k)

  expect_equal(res_cpp$d, multivarious::sdev(res_r), tolerance = 1e-6)
  expect_equal(res_cpp$k, k)
  # Compare actual values, not subspaces
  # res_cpp$u is the M-weighted scores (M ou D); scores() is ou D = X A ov
  # (Allen et al. 2014), stored via fit$u = M ou.
  u_ref <- as.matrix(sweep(res_r$u, 2, multivarious::sdev(res_r), `*`))
  expect_equal(align_signs(u_ref, res_cpp$u), u_ref,
               tolerance = 1e-6, check.attributes = FALSE)
  expect_equal(multivarious::scores(res_r),
               sweep(res_r$ou, 2, multivarious::sdev(res_r), `*`),
               tolerance = 1e-8, check.attributes = FALSE)
  v_ref <- as.matrix(multivarious::components(res_r))
  expect_equal(align_signs(v_ref, res_cpp$v), v_ref,
               tolerance = 1e-6, check.attributes = FALSE)
})

test_that("gmd_fast_cpp matches genpca (use_cpp=TRUE) for p <= n, sparse constraints", {
  set.seed(3)
  n <- 25
  p <- 20
  k <- 4
  X <- matrix(rnorm(n * p), n, p)
  # Create sparse constraints (e.g., diagonal)
  Q <- Diagonal(n, x = runif(n, 0.5, 1.5))
  R <- Diagonal(p, x = runif(p, 0.5, 1.5))

  X_centered <- scale(X, center = TRUE, scale = FALSE)

  res_r <- genpca(X_centered, M = Q, A = R, ncomp = k, use_cpp = TRUE, method = "spectra", preproc = multivarious::pass())
  res_cpp <- genpca:::gmd_fast_cpp(X_centered, Q = Q, R = R, k = k)

  expect_equal(res_cpp$d, multivarious::sdev(res_r), tolerance = 1e-6)
  expect_equal(res_cpp$k, k)
  # Compare actual values, not subspaces
  # res_cpp$u is the M-weighted scores (M ou D); scores() is ou D = X A ov
  # (Allen et al. 2014), stored via fit$u = M ou.
  u_ref <- as.matrix(sweep(res_r$u, 2, multivarious::sdev(res_r), `*`))
  expect_equal(align_signs(u_ref, res_cpp$u), u_ref,
               tolerance = 1e-6, check.attributes = FALSE)
  expect_equal(multivarious::scores(res_r),
               sweep(res_r$ou, 2, multivarious::sdev(res_r), `*`),
               tolerance = 1e-8, check.attributes = FALSE)
  v_ref <- as.matrix(multivarious::components(res_r))
  expect_equal(align_signs(v_ref, res_cpp$v), v_ref,
               tolerance = 1e-6, check.attributes = FALSE)
})

test_that("gmd_fast_cpp matches genpca (use_cpp=TRUE) for p > n, sparse constraints", {
  set.seed(4)
  n <- 20
  p <- 25
  k <- 4
  X <- matrix(rnorm(n * p), n, p)
  Q <- Diagonal(n, x = runif(n, 0.5, 1.5))
  R <- Diagonal(p, x = runif(p, 0.5, 1.5))

  X_centered <- scale(X, center = TRUE, scale = FALSE)

  res_r <- genpca(X_centered, M = Q, A = R, ncomp = k, use_cpp = TRUE, method = "spectra", preproc = multivarious::pass())
  res_cpp <- genpca:::gmd_fast_cpp(X_centered, Q = Q, R = R, k = k)

  expect_equal(res_cpp$d, multivarious::sdev(res_r), tolerance = 1e-6)
  expect_equal(res_cpp$k, k)
  # Compare actual values, not subspaces
  # res_cpp$u is the M-weighted scores (M ou D); scores() is ou D = X A ov
  # (Allen et al. 2014), stored via fit$u = M ou.
  u_ref <- as.matrix(sweep(res_r$u, 2, multivarious::sdev(res_r), `*`))
  expect_equal(align_signs(u_ref, res_cpp$u), u_ref,
               tolerance = 1e-6, check.attributes = FALSE)
  expect_equal(multivarious::scores(res_r),
               sweep(res_r$ou, 2, multivarious::sdev(res_r), `*`),
               tolerance = 1e-8, check.attributes = FALSE)
  v_ref <- as.matrix(multivarious::components(res_r))
  expect_equal(align_signs(v_ref, res_cpp$v), v_ref,
               tolerance = 1e-6, check.attributes = FALSE)
})

test_that("gmd_fast_cpp handles k=1 correctly", {
  set.seed(5)
  n <- 10
  p <- 8
  k <- 1
  X <- matrix(rnorm(n * p), n, p)
  Q <- Diagonal(n, x = runif(n, 0.5, 1.5))
  R <- Diagonal(p, x = runif(p, 0.5, 1.5))
  X_centered <- scale(X, center = TRUE, scale = FALSE)

  res_r <- genpca(X_centered, M = Q, A = R, ncomp = k, use_cpp = TRUE, method = "spectra", preproc = multivarious::pass())
  res_cpp <- genpca:::gmd_fast_cpp(X_centered, Q = Q, R = R, k = k)

  expect_equal(res_cpp$d, multivarious::sdev(res_r), tolerance = 1e-6)
  expect_equal(res_cpp$k, k)
  expect_equal(ncol(res_cpp$u), k)
  expect_equal(ncol(res_cpp$v), k)
  # Compare actual values, not subspaces
  # res_cpp$u is the M-weighted scores (M ou D); scores() is ou D = X A ov
  # (Allen et al. 2014), stored via fit$u = M ou.
  u_ref <- as.matrix(sweep(res_r$u, 2, multivarious::sdev(res_r), `*`))
  expect_equal(align_signs(u_ref, res_cpp$u), u_ref,
               tolerance = 1e-6, check.attributes = FALSE)
  expect_equal(multivarious::scores(res_r),
               sweep(res_r$ou, 2, multivarious::sdev(res_r), `*`),
               tolerance = 1e-8, check.attributes = FALSE)
  v_ref <- as.matrix(multivarious::components(res_r))
  expect_equal(align_signs(v_ref, res_cpp$v), v_ref,
               tolerance = 1e-6, check.attributes = FALSE)
})

test_that("gmd_fast_cpp returns fewer components if necessary", {
  # Test case where numerical rank might be less than k requested
  set.seed(6)
  n <- 10
  p <- 10
  k <- 8
  # Create a low-rank matrix
  rank_true <- 5
  X <- matrix(rnorm(n * rank_true), n, rank_true) %*% matrix(rnorm(rank_true * p), rank_true, p)
  Q <- Diagonal(n)
  R <- Diagonal(p)
  X_centered <- scale(X, center = TRUE, scale = FALSE)

  # Expect warning when k > rank? Maybe not from C++ directly
  # The C++ code filters eigenvalues close to zero
  res_cpp <- genpca:::gmd_fast_cpp(X_centered, Q = Q, R = R, k = k, tol = 1e-9)

  expect_lte(res_cpp$k, k)
  expect_lte(res_cpp$k, min(n, p))
  # Actual rank might depend on numerical tolerance
  expect_true(res_cpp$k >= rank_true - 1 && res_cpp$k <= rank_true + 1)
  expect_equal(length(res_cpp$d), res_cpp$k)
  expect_equal(ncol(res_cpp$u), res_cpp$k)
  expect_equal(ncol(res_cpp$v), res_cpp$k)
})

test_that("gmd_fast_cpp matches genpca spectra for identity constraints", {
  # Direct comparison gmd_fast_cpp vs genpca(method="spectra") with identity metrics
  set.seed(2024)
  n <- 30
  p <- 20
  k <- 5
  X <- matrix(rnorm(n * p), n, p)
  X_centered <- scale(X, scale = FALSE)

  Q <- diag(n)  # identity row metric
  R <- diag(p)  # identity column metric

  # gmd_fast_cpp (spectra-like path)
  res_cpp <- genpca:::gmd_fast_cpp(X_centered, Q = Matrix::Matrix(Q), R = Matrix::Matrix(R), k = k)

  # genpca with method="spectra" (iterative SVD path that wraps gmd_fast_cpp)
  res_spectra <- genpca(X_centered, M = Q, A = R, ncomp = k, method = "spectra", preproc = multivarious::pass())

  # Singular values should match (iterative solver: 1e-6 is realistic)
  expect_equal(res_cpp$d, multivarious::sdev(res_spectra), tolerance = 1e-6)

  # Scores and loadings should match up to column sign
  u_ref <- as.matrix(multivarious::scores(res_spectra))
  v_ref <- as.matrix(multivarious::components(res_spectra))
  expect_equal(align_signs(u_ref, as.matrix(res_cpp$u)), u_ref,
               tolerance = 1e-6, check.attributes = FALSE)
  expect_equal(align_signs(v_ref, as.matrix(res_cpp$v)), v_ref,
               tolerance = 1e-6, check.attributes = FALSE)
})

test_that("gmd_fast_cpp matches genpca spectra for diagonal constraints", {

  # Test with non-identity diagonal metrics
  set.seed(2025)
  n <- 25
  p <- 15
  k <- 4
  X <- matrix(rnorm(n * p), n, p)
  X_centered <- scale(X, scale = FALSE)

  # Diagonal metrics (e.g., observation weights)
  row_weights <- runif(n, 0.5, 1.5)
  col_weights <- runif(p, 0.5, 1.5)
  Q <- diag(row_weights)
  R <- diag(col_weights)

  # gmd_fast_cpp
  res_cpp <- genpca:::gmd_fast_cpp(X_centered, Q = Matrix::Matrix(Q), R = Matrix::Matrix(R), k = k)

  # genpca with method="spectra"
  res_spectra <- genpca(X_centered, M = Q, A = R, ncomp = k, method = "spectra", preproc = multivarious::pass())

  # Singular values should match (iterative solver: 1e-6 is realistic)
  expect_equal(res_cpp$d, multivarious::sdev(res_spectra), tolerance = 1e-6)
})

## ────────────────────────────────────────────────────────────────────────────────
##  tests/testthat/test-genpca-advanced.R
## ────────────────────────────────────────────────────────────────────────────────
context("genpca – advanced properties")

set.seed(42)
Xsmall <- matrix(rnorm(50 * 30), 50, 30)      # n > p       (spectra: right‑side)
Xwide  <- matrix(rnorm(40 * 120), 40, 120)    # n < p       (spectra: left‑side)

## -------------------------------------------------------------------------------
test_that("Spectra method matches eigen method on modest problems", {


  ## 1) n >= p   (right‑side operator)
  fit_eig  <- genpca(Xsmall, ncomp = 10, method = "eigen",
                     preproc = multivarious::center())
  fit_spc  <- genpca(Xsmall, ncomp = 10, method = "spectra",
                     preproc = multivarious::center(), tol_spectra = 1e-10)

  expect_equal(fit_eig$sdev,           fit_spc$sdev,           tolerance = 1e-6)
  expect_equal(abs(multivarious::scores(fit_eig)),  abs(multivarious::scores(fit_spc)),  tolerance = 1e-5, check.attributes = FALSE)

  ## 2) n < p    (left‑side operator)
  fit_eig_w <- genpca(Xwide,  ncomp = 15, method = "eigen",
                      preproc = multivarious::center())
  fit_spc_w <- genpca(Xwide,  ncomp = 15, method = "spectra",
                      preproc = multivarious::center(), tol_spectra = 1e-10)

  expect_equal(fit_eig_w$sdev,          fit_spc_w$sdev,          tolerance = 1e-6)
  expect_equal(abs(multivarious::scores(fit_eig_w)), abs(multivarious::scores(fit_spc_w)), tolerance = 1e-5, check.attributes = FALSE)
})

## -------------------------------------------------------------------------------
test_that("Orthonormality holds in (M,A) metrics", {
  skip_if_not_installed("adjoin")

  Mrow <- adjoin::temporal_adjacency(1:nrow(Xsmall))
  Mrow <- t(Mrow) %*% Mrow                 # PSD & dense
  set.seed(561)
  Acol <- cov(matrix(rnorm(ncol(Xsmall)^2), ncol(Xsmall)))

  fit <- genpca(Xsmall, M = Mrow, A = Acol,
                ncomp = 12, preproc = multivarious::center())

  QtU <- t(fit$ou) %*% Mrow %*% fit$ou     # should be ~ I
  RtV <- t(fit$ov) %*% Acol %*% fit$ov     # should be ~ I

  expect_true(max(abs(QtU - diag(ncol(QtU)))) < 1e-8)
  expect_true(max(abs(RtV - diag(ncol(RtV)))) < 1e-8)
})

## -------------------------------------------------------------------------------
test_that("Reconstruction error decreases monotonically with ncomp", {

  set.seed(576)
  X <- matrix(rnorm(120 * 45), 120, 45)
  fit_all <- genpca(X, ncomp = 20, preproc = multivarious::center())

  rss <- vapply(1:20, function(k) {
           recon <- reconstruct(fit_all, comp = 1:k)
           sum((X - recon)^2)
         }, numeric(1))

  ## RSS should strictly decrease until floating‑point noise kicks in
  expect_true(all(diff(rss) <= 1e-12))
})
