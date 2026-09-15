context("gmd_deflationR empty case")

library(Matrix)
library(multivarious)

# gmd_deflationR is internal

test_that("gmd_deflationR returns zero-sized matrices when no components found", {
  X <- matrix(0, 5, 4)
  Q <- Diagonal(5)
  R <- Diagonal(4)
  res <- genpca:::gmd_deflationR(X, Q, R, k = 3, verbose = FALSE)
  expect_equal(res$k, 0)
  expect_equal(dim(res$u), c(5, 0))
  expect_equal(dim(res$v), c(4, 0))
  expect_equal(length(res$d), 0)
  expect_equal(res$u, matrix(0, 5, 0))
  expect_equal(res$v, matrix(0, 4, 0))
})


test_that("genpca handles zero-rank input consistently", {
  X <- matrix(0, 6, 4)
  fit <- genpca(X, ncomp = 3, method = "deflation", use_cpp = FALSE)
  expect_s3_class(fit, c("genpca", "bi_projector", "projector"))
  expect_equal(multivarious::ncomp(fit), 0)
  expect_equal(dim(multivarious::scores(fit)), c(6, 0))
  expect_equal(dim(fit$v), c(4, 0))
  expect_equal(length(multivarious::sdev(fit)), 0)
})
