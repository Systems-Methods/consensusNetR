test_that("get_gene_community_membership works", {
  expect_snapshot(
    get_gene_community_membership(testing_consensus_comms, testing_memb_list)
  )
})
