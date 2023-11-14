test_that("plot_consensus_eig_dist defaults", {
  save_png <- function(code, width = 300, height = 300) {
    path <- withr::local_tempfile(fileext = ".png")
    png(path, width = width, height = height, pointsize = 8)
    on.exit(dev.off())
    code

    path
  }

  test_consensus_memb <- suppressMessages(calc_consensus_memberships(
    testing_consensus_comms,
    testing_memb_list,
    gene_cohort_N = 2
  ))

  eigen_list <- suppressMessages(list(
    GSE39582_eigen = icWGCNA::compute_eigengene_matrix(
      testing_ex_list$GSE39582, test_consensus_memb
    ),
    READ_eigen = icWGCNA::compute_eigengene_matrix(
      testing_ex_list$READ, test_consensus_memb
    ),
    COAD_eigen = icWGCNA::compute_eigengene_matrix(
      testing_ex_list$COAD, test_consensus_memb
    )
  ))

  # only tested exact figure in Mac OS
  if (Sys.info()["sysname"] == "Darwin") {
    expect_snapshot_file(
      save_png(plot_consensus_eig_dist(eigen_list)),
      "plot_consensus_eig_dist.png"
    )
  } else if (Sys.info()["sysname"] == "Linux") {
    expect_no_error(
      plot_consensus_eig_dist(eigen_list)
    )
  }
  expect_no_error(
    plot_consensus_eig_dist(
      eigen_list,
      filename = withr::local_tempfile(fileext = ".png"), dpi = 10
    )
  )
})
