context("Number of stable generalized eigenvalues")

test_that("num_stable with beta = NULL", {
  expect_equal(num_stable(c(1)), 1)
  expect_equal(num_stable(c(1.1 + 0i)), 0)
  expect_equal(num_stable(c(0.9, 1.1i, 0.8)), 2)
})

test_that("num_stable with beta != NULL", {
  expect_equal(num_stable(c(1), c(2)), 1)
  expect_equal(num_stable(c(1.1 + 0i), c(0.8)), 0)
  expect_equal(num_stable(c(0.9, 1.1i, 0.8), c(0.9, 1.1i, 0.8)), 3)
})
