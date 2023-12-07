#' Compute RBH (reciprocal best hits) for Two Networks based on Overlap of top genes
#'
#' Compute RBH (reciprocal best hits) for Two Networks based on Overlap of top
#' community genes. They must have row names containing some sort of common geneID
#' vocabulary, so that they can be intersected to get common rows.
#'
#' @param net_membership_1 community membership matrix from first network.
#' @param net_membership_2 community membership matrix from second network.
#' @inheritParams construct_rbh_overlap_based
#'
#' @return Matrix with rows matching the columns in net_membership_1 and columns
#' matching the columns of net_membership_2, with values 1 indicating two metagenes are a
#' hit, and 0 indicating that they are not.
#' @keywords internal
#' @examples
#' \dontrun{
#' rbh <- construct_2study_rbh_overlap_based(
#'   GSE39582_icwgcna$community_membership,
#'   read_icwgcna$community_membership
#' )
#' }
construct_2study_rbh_overlap_based <- function(net_membership_1,
                                               net_membership_2,
                                               top_n = 50,
                                               memb_cut = 0) {
  ret <- matrix(0, ncol(net_membership_1), ncol(net_membership_2))
  net_membership_1 <- net_membership_1[rownames(net_membership_1) %in%
    rownames(net_membership_2), ]
  net_membership_2 <- net_membership_2[match(
    rownames(net_membership_1),
    rownames(net_membership_2)
  ), ]
  rankX <- apply(
    apply(-net_membership_1, 2, rank) <= top_n &
      net_membership_1 >= memb_cut,
    2, as.numeric
  )
  rankY <- apply(
    apply(-net_membership_2, 2, rank) <= top_n &
      net_membership_2 >= memb_cut,
    2, as.numeric
  )

  overlap <- t(rankX) %*% rankY

  x2y <- max.col(overlap, ties.method = "first")
  x2y[rowSums(overlap) == 0] <- NA
  x2y <- data.frame(
    x2yName = 1:length(x2y),
    x2yMap = x2y
  )
  y2x <- max.col(t(overlap), ties.method = "first")
  y2x[rowSums(t(overlap)) == 0] <- NA
  y2x <- data.frame(
    y2xName = 1:length(y2x),
    y2xMap = y2x
  )
  y2x <- y2x[match(x2y$x2yMap, y2x$y2xName), ]

  edges <- as.matrix(x2y[which(x2y$x2yName == y2x$y2xMap), ])
  ret[edges] <- overlap[edges]
  ret / top_n
}


#' Compute correlation based RBH for Two Metagene Datasets
#'
#' Computes a single network with two datasets. They must have rownames
#' containing some sort of common geneID vocabulary, so that they can
#' be intersected to get common rows.
#'
#' If one wants the reciprocal best hits network to have separate dimension
#' names for rows and columns, the column names for net_membership_1 and net_membership_2 should
#' be distinct.
#'
#' @param net_membership_1 community membership matrix from first network.
#' @param net_membership_2 community membership matrix from second network.
#' @inheritParams construct_rbh_correlation_based
#'
#' @return Matrix with rows matching the columns in net_membership_1 and columns
#' matching the columns of net_membership_2, with values 1 indicating two
#' metagenes are a hit, and 0 indicating that they are not.
#' @keywords internal
#' @examples
#' \dontrun{
#' rbh <- construct_2study_rbh_correlation_based(
#'   GSE39582_icwgcna$community_membership,
#'   read_icwgcna$community_membership
#' )
#' }
construct_2study_rbh_correlation_based <- function(
    net_membership_1,
    net_membership_2,
    lower_quant = 0,
    upper_quant = 1.0,
    max_rank = 1,
    abs = FALSE,
    sparse = FALSE,
    method = "pearson") {
  # Filter for rows with common gene names
  net_membership_1 <- stats::na.omit(net_membership_1)
  net_membership_2 <- stats::na.omit(net_membership_2)
  common_rows <- intersect(
    row.names(net_membership_1),
    row.names(net_membership_2)
  )
  net_membership_1 <- net_membership_1[common_rows, ]
  net_membership_2 <- net_membership_2[common_rows, ]

  # Calculate correlation matrix for the columns
  if (abs) {
    rbh <- abs(stats::cor(net_membership_1, net_membership_2, method = method))
  } else {
    rbh <- stats::cor(net_membership_1, net_membership_2, method = method)
  }
  gen_cutoff <- stats::quantile(rbh, upper_quant, na.rm = TRUE)
  rbh_threshold <- stats::quantile(rbh, lower_quant, na.rm = TRUE)


  # Construct a list representing the best hit(s) for each row (up to max_rank
  # hits), and a list representing the same for each column.
  row_count <- nrow(rbh)
  col_count <- ncol(rbh)
  if (max_rank == 1) {
    # Use relatively fast argmax (max_col)
    row_hits <- max.col(rbh)
    col_hits <- max.col(t(rbh))
  } else {
    # Use slower argsort function (order)
    rows <- lapply(seq_len(row_count), function(i) rbh[i, ])
    cols <- lapply(seq_len(col_count), function(i) rbh[, i])
    row_hits <- lapply(rows, function(row) {
      return(order(row)[(row_count - max_rank + 1):row_count])
    })
    col_hits <- lapply(cols, function(col) {
      return(order(col)[(col_count - max_rank + 1):col_count])
    })
  }

  # Construct a list of lists containing booleans representing whether
  # (net_membership_1[, i], net_membership_2[, j]) is a reciprocal best hit, based on our
  # parameters.
  temp <- lapply(1:row_count, function(row_num) {
    lapply(1:col_count, function(col_num) {
      return(ifelse((rbh[row_num, col_num] >= gen_cutoff) |
        (col_num %in% row_hits[[row_num]] &
          row_num %in% col_hits[[col_num]] &
          rbh[row_num, col_num] >= rbh_threshold),
      rbh[row_num, col_num],
      0
      ))
    })
  })
  # Bind these lists into a matrix
  rbh_network <- Matrix::Matrix(as.numeric(do.call(rbind, temp)),
    nrow = row_count,
    sparse = sparse,
    dimnames = list(
      dimnames(net_membership_1)[[2]],
      dimnames(net_membership_2)[[2]]
    )
  )
  return(as.matrix(rbh_network))
}
