test_that("calc_cor_of_cor works", {
  expect_snapshot(calc_cor_of_cor(testing_ex_list))
})

test_that("calc_cor_of_cor parallel", {
  local_mocked_bindings(getDoParName = function(package) "doMC")
  expect_message(
    suppressWarnings(calc_cor_of_cor(testing_ex_list)),
    paste0("Using doMC with ", foreach::getDoParWorkers(), " workers"))
})


test_that("compare_networks works", {
  expect_snapshot(suppressWarnings(compare_networks(testing_memb_list[[1]],
                                   testing_memb_list[[2]])))
})

