test_that("get_gene_community_membership works", {
  test_rbh <- suppressMessages(construct_rbh_overlap_based(testing_memb_list))
  test_consensus <- suppressMessages(detect_consensus_communities(test_rbh))

  expect_snapshot(
    get_gene_community_membership(test_consensus, testing_memb_list)
  )
})
