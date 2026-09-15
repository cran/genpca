library(testthat)

# Additional tests for rpls

test_that("rpls errors when lambda vector shorter than K", {
  set.seed(1)
  X <- matrix(rnorm(10 * 3), 10, 3)
  Y <- matrix(rnorm(10 * 2), 10, 2)
  lambda_vec <- c(0.1, 0.2)  # length 2, but K=3
  expect_error(rpls(X, Y, K = 3, lambda = lambda_vec), "lambda must be either")
})

test_that("nonneg flag warns when penalty is ridge", {
  set.seed(2)
  X <- matrix(rnorm(8 * 4), 8, 4)
  Y <- matrix(rnorm(8 * 2), 8, 2)
  expect_warning(rpls(X, Y, K = 1, penalty = "ridge", nonneg = TRUE),
                 "nonneg option is ignored")
})
