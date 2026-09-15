test_that("repair_metric leaves PSD input unchanged and reports it", {
  A <- crossprod(matrix(rnorm(30), 10, 3))          # PSD, rank 3
  B <- repair_metric(A, method = "ridge")
  rep <- attr(B, "repair_report")
  expect_s3_class(rep, "metric_repair_report")
  expect_false(rep$changed)
  expect_equal(rep$shift, 0)
  expect_equal(rep$rank, 3L)
  expect_equal(as.matrix(B), A, tolerance = 1e-12, ignore_attr = TRUE)
  expect_output(print(rep), "changed:")
})

test_that("repair_metric ridge and clip both return PSD matrices with a diagnostic", {
  A <- matrix(c(1, 2, 2, 1), 2)                    # eigenvalues 3, -1
  Br <- repair_metric(A, method = "ridge")
  rr <- attr(Br, "repair_report")
  expect_true(rr$changed)
  expect_equal(rr$min_eigenvalue_before, -1)
  expect_true(rr$min_eigenvalue_after > 0)
  expect_true(is.finite(rr$shift) && rr$shift > 1)
  expect_equal(as.matrix(Br), A + diag(rr$shift, 2), tolerance = 1e-12, ignore_attr = TRUE)

  Bc <- repair_metric(A, method = "clip")
  rc <- attr(Bc, "repair_report")
  expect_true(rc$changed)
  expect_true(is.na(rc$shift))
  ev <- eigen(as.matrix(Bc), symmetric = TRUE)$values
  expect_equal(ev, c(3, 0), tolerance = 1e-12)
  expect_equal(rc$rank, 1L)
})

test_that("repair_metric rejects asymmetric input and handles sparse input", {
  expect_error(repair_metric(matrix(c(1, 0.5, 0, 1), 2)), "symmetric")
  n <- 50
  S <- Matrix::bandSparse(n, n, k = c(-1, 0, 1),
                          diagonals = list(rep(-1, n - 1), rep(1, n), rep(-1, n - 1)),
                          symmetric = FALSE)
  S <- Matrix::forceSymmetric(S)                    # indefinite tridiagonal
  B <- repair_metric(S, method = "ridge", name = "M")
  rep <- attr(B, "repair_report")
  expect_true(rep$changed)
  expect_equal(rep$name, "M")
  expect_true(methods::is(B, "sparseMatrix"))
  expect_true(genpca:::is_pd(B))
})
