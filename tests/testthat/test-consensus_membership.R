test_that("calc_consensus_memberships works", {
  expect_snapshot(calc_consensus_memberships(testing_consensus_comms,
    testing_memb_list,
    gene_cohort_N = 2
  ))
})

test_that("calc_consensus_memberships more param", {
  expect_snapshot(calc_consensus_memberships(testing_consensus_comms,
    testing_memb_list,
    compressIntra = "mean",
    weights = "Dynamic"
  ))
})

test_that("calc_consensus_memberships weights", {
  expect_snapshot(calc_consensus_memberships(testing_consensus_comms,
    testing_memb_list,
    weights = c(.2, .2, .6)
  ))
})

test_that("calc_consensus_memberships weights", {
  expect_error(
    calc_consensus_memberships(testing_consensus_comms,
                               testing_memb_list,
                               gene_cohort_N = 10),
    "gene_cohort_N can't be larger than number of studies in network_membership_list"
  )
})

