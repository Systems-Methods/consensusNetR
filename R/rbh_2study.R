#' Compute RBH (reciprocal best hits) for Two Networks based on Overlap of top genes
#'
#' Compute RBH (reciprocal best hits) for Two Networks based on Overlap of top
#' community genes. They must have row names containing some sort of common geneID
#' vocabulary, so that they can be intersected to get common rows.
#'
#' @param meta_g1 community membership matrix from first network.
#' @param meta_g2 community membership matrix from second network.
#' Row names of these matrices must be gene names, and the column names should be
#' unique community names.
#' @param top_n number of top genes to usin in compuitng overlap
#' @param memb_cut membership threshold for stricter thresholding. Only genes with
#' membership score greater the threshold are used
#'
#' @return Matrix with rows matching the columns in meta_g1 and columns
#' matching the columns of meta_g2, with values 1 indicating two metagenes are a
#' hit, and 0 indicating that they are not.
#'
#'
construct_2study_rbh_overlap_based <- function(meta_g1,
                                               meta_g2,
                                               top_n = 50,
                                               memb_cut = 0) {
  #simpler approach that doesn't punish the overlap measure for genes not shared by platforms
  ret        <- matrix(0, ncol(meta_g1), ncol(meta_g2))
  meta_g1      <- meta_g1[rownames(meta_g1) %in% rownames(meta_g2),]
  meta_g2      <- meta_g2[match(rownames(meta_g1),rownames(meta_g2)),]
  rankX      <- apply(apply(-meta_g1,2,rank) <= top_n & meta_g1 >= memb_cut,
                      2, as.numeric)
  rankY      <- apply(apply(-meta_g2,2,rank) <= top_n & meta_g2 >= memb_cut,
                      2, as.numeric)

  overlap    <- t(rankX) %*% rankY
  rm(meta_g1,meta_g2, rankX,rankY)
  x2y        <- plyr::alply(overlap, 1,
                            function(x) {
                              ret <- which(x == max(x) & x != 0);
                              if (length(ret) == 0) {
                                return(NA)
                              } else {
                                return(ret[1])
                              }});
  y2x        <- plyr::alply(overlap, 2,
                            function(x) {
                              ret <- which(x == max(x) & x != 0)
                              if (length(ret) == 0) {
                                return(NA)
                              } else {
                                return(ret[1])
                              }});
  x2y        <- data.frame(x2yName = as.numeric(names(x2y)),
                           x2yMap = unlist(x2y))
  y2x        <- data.frame(y2xName = as.numeric(names(y2x)),
                           y2xMap = unlist(y2x))
  y2x        <- y2x[match(x2y$x2yMap,y2x$y2xName),]

  edges      <- as.matrix(x2y[which(x2y$x2yName == y2x$y2xMap),])
  ret[edges] <- overlap[edges]
  ret        <- ret/top_n
  return(ret)
}


#' Compute correlation based RBH for Two Metagene Datasets
#'
#' Computes a single network with two datasets. They must have rownames
#' containing some sort of common geneID vocabulary, so that they can
#' be intersected to get common rows.
#'
#' If one wants the reciprocal best hits network to have separate dimension
#' names for rows and columns, the column names for meta_g1 and meta_g2 should
#' be distinct.
#'
#' @param meta_g1,meta_g2 metagene matrices from two different datasets. Row
#' names of these matrices must be gene names, and the column names should be
#' unique metagene names.
#' @inheritParams construct_rbh_correlation_based
#'
#' @return Matrix with rows matching the columns in meta_g1 and columns
#' matching the columns of meta_g2, with values 1 indicating two metagenes are a
#' hit, and 0 indicating that they are not.
#'
#' @examples
#' \dontrun{
#' network_file_dir <- system.file("extdata", package = "consensusNetR")
#' network_file_list <- list.files(network_file_dir, full.names = TRUE)
#' ds1 <- as.matrix(data.table::fread(network_file_list[[1]]))
#' row.names(ds1) <- ds1[, 1]
#' ds1 <- ds1[, -1]
#' colnames(ds1) <- paste0(colnames(ds1), sub(
#'   ".*/", "",
#'   network_file_list[[1]]
#' ))
#'
#' ds2 <- as.matrix(data.table::fread(network_file_list[[2]]))
#' row.names(ds2) <- ds2[, 1]
#' ds2 <- ds2[, -1]
#' colnames(ds2) <- paste0(colnames(ds2), sub(
#'   ".*/", "",
#'   network_file_list[[2]]
#' ))
#' rbh <- construct_2study_rbh_correlation_based(ds1, ds2)
#' }
#'
construct_2study_rbh_correlation_based <- function(meta_g1, meta_g2, lower_quant = 0,
                                                 upper_quant = 1.0, max_rank = 1,
                                                 abs = FALSE, sparse = FALSE, method = "pearson") {

  # Filter for rows with common gene names
  meta_g1 <- stats::na.omit(meta_g1)
  meta_g2 <- stats::na.omit(meta_g2)
  common_rows <- intersect(row.names(meta_g1), row.names(meta_g2))
  meta_g1 <- meta_g1[common_rows, ]
  meta_g2 <- meta_g2[common_rows, ]

  # Calculate correlation matrix for the columns
  if (abs) {
    rbh <- abs(stats::cor(meta_g1, meta_g2, method = method))
  } else {
    rbh <- stats::cor(meta_g1, meta_g2, method = method)
  }
  gen_cutoff <- stats::quantile(rbh, upper_quant,na.rm = TRUE)
  rbh_threshold <- stats::quantile(rbh, lower_quant,na.rm = TRUE)


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
  # (meta_g1[, i], meta_g2[, j]) is a reciprocal best hit, based on our
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
                                  dimnames(meta_g1)[[2]],
                                  dimnames(meta_g2)[[2]]
                                )
  )
  return(as.matrix(rbh_network))
}
