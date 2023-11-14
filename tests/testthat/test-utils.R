test_that("fisherZ check", {
  expect_equal(fisherZ(c(0, 1)), c(0, Inf))
  expect_equal(fisherZ(c(.3, .7)), c(0.30951960, 0.86730053))
})
