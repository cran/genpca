library(testthat)
library(multivarious)
library(Matrix)

context("reconstruct with zero components")

test_that("reconstruct returns matrix with expected dimensions when ncomp is zero", {
  n <- 5
  p <- 4
  A <- Matrix::Diagonal(p)
  M <- Matrix::Diagonal(n)
  # Create a proper pre_processor
  set.seed(13)
  dummy_data <- matrix(rnorm(n * p), n, p)
  dummy_fit <- genpca(dummy_data, preproc = multivarious::pass())
  proper_preproc <- dummy_fit$preproc

  obj <- multivarious::bi_projector(
    v = matrix(0, p, 0),
    s = matrix(0, n, 0),
    sdev = numeric(0),
    preproc = proper_preproc,
    ov = matrix(0, p, 0),
    ou = matrix(0, n, 0),
    u  = matrix(0, n, 0),
    classes = "genpca",
    A = A,
    M = M,
    propv = numeric(0),
    cumv = numeric(0)
  )

  recon_default <- reconstruct(obj)
  expect_equal(dim(recon_default), c(n, p))

  ri <- 2:4
  ci <- c(1, 3)
  recon_subset <- reconstruct(obj, rowind = ri, colind = ci)
  expect_equal(dim(recon_subset), c(length(ri), length(ci)))
})
