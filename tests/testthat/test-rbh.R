test_that("rbh overlap defaults", {
  expect_snapshot(construct_rbh_overlap_based(testing_memb_list))
})
test_that("rbh overlap other params", {
  expect_snapshot(construct_rbh_overlap_based(testing_memb_list,
    memb_cut = .9
  ))
})

test_that("rbh correlation defaults", {
  expect_snapshot(construct_rbh_correlation_based(testing_memb_list))
})
test_that("rbh correlation other params", {
  expect_snapshot(construct_rbh_correlation_based(testing_memb_list,
    abs = TRUE,
    max_rank = 2
  ))
})

test_that("plot rbh defaults", {
  # only tested exact figure in Mac OS
  if (Sys.info()["sysname"] == "Darwin") {
    save_png <- function(code, width = 400, height = 400) {
      path <- withr::local_tempfile(fileext = ".png")
      png(path, width = width, height = height)
      on.exit(dev.off())
      code

      path
    }

    test_rbh <- suppressMessages(construct_rbh_overlap_based(testing_memb_list))
    expect_snapshot_file(
      save_png(plot_rbh(test_rbh, testing_memb_list)),
      "plot_rbh.png"
    )
  }
})


test_that("plot rbh linux", {
  # only tested exact figure in Mac OS
  if (Sys.info()["sysname"] == "Linux") {
    save_png <- function(code, width = 400, height = 400) {
      path <- withr::local_tempfile(fileext = ".png")
      png(path, width = width, height = height)
      on.exit(dev.off())
      code

      path
    }

    test_rbh <- suppressMessages(construct_rbh_overlap_based(testing_memb_list))
    expect_snapshot_file(
      save_png(plot_rbh(test_rbh, testing_memb_list)),
      "plot_rbh_linux.png"
    )
  }
})

