context("QZ decomposition")

test_that("Hansen RBC", {
  A <- structure(
    c(0, 0.33, 1.01010101010101, 0, 0, 1, 0.100306091215182,
      0.95, 0, -0.82676426385235, 0.067205081114172, 0,
      -1.01010101010101, -1, -0.0773060912151821, 0),
    .Dim = c(4L, 4L))
  E <- structure(
    c(-0.0221776767676768, 0, 1, 0, 0.0331010101010101,
      0, 0, 1, 0.0221776767676768, 0, 0, 0, -1.01010101010101,
      0, 0, 0), .Dim = c(4L, 4L))
  eig <- c(0.951445532228811, 0.95, 1.0616488026748, Inf)
  ret <- qz(A, E)

  testthat::expect_equal(sort(ret$ALPHA / ret$BETA), sort(eig))
  testthat::expect_equal(ret$ALPHA / ret$BETA <= 1,
                         c(TRUE, TRUE, FALSE, FALSE))
  testthat::expect_equal(t(ret$Q) %*% A %*% ret$Z, ret$T)
  testthat::expect_equal(t(ret$Q) %*% E %*% ret$Z, ret$S)
})
