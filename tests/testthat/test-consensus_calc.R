test_that("detect_consensus_communities works", {
  test_rbh <- suppressMessages(construct_rbh_overlap_based(testing_memb_list))
  expect_snapshot(detect_consensus_communities(test_rbh))
})
