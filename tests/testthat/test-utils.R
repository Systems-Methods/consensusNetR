test_that("fisherZ check", {
  expect_equal(fisherZ(c(0, 1)), c(0, Inf))
  expect_equal(fisherZ(c(.3, .7)), c(0.30951960, 0.86730053))
})

test_that("check_installed", {
  local_mocked_bindings(is_installed = function(package) FALSE)
  expect_error(
    check_installed('weird_pkg'),
    "Must have the following R packages installed for this function: weird_pkg")
  expect_error(
    check_installed('stats'),
    "Must have the following R packages installed for this function: stats")
})

