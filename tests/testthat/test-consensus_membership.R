test_that("calc_consensus_memberships works", {
  test_rbh <- suppressMessages(construct_rbh_overlap_based(testing_memb_list))
  test_consensus <- suppressMessages(detect_consensus_communities(test_rbh))
  expect_snapshot(calc_consensus_memberships(test_consensus, testing_memb_list,
                                             gene_cohort_N = 2))
})

test_that("calc_consensus_memberships more param", {
  test_rbh <- suppressMessages(construct_rbh_overlap_based(testing_memb_list))
  test_consensus <- suppressMessages(detect_consensus_communities(test_rbh))
  expect_snapshot(calc_consensus_memberships(test_consensus, testing_memb_list,
                                             compressIntra = "mean",
                                             weights = "Dynamic"))
})

test_that("calc_consensus_memberships weights", {
  test_rbh <- suppressMessages(construct_rbh_overlap_based(testing_memb_list))
  test_consensus <- suppressMessages(detect_consensus_communities(test_rbh))
  expect_snapshot(calc_consensus_memberships(test_consensus, testing_memb_list,
                                             weights = c(.2,.2,.6)))
})
