#' Construct Meta Reciprocal Best Hits based on overlaps
#'
#' Iterate through a list containing membership matrices and construct pairwise
#' networks for every combination, and output them as a symmetric matrix with
#' labeled dimensions. This can be interpreted as an "adjacency" matrix between
#' all of the different metagenes across networks in the input list. The
#' column/row names will take each networks communities' names and concatenate them
#' with the network list names, such that the dimension names are guaranteed to
#' be unique.
#'
#' @param network_membership_list a list containing community membership scores for each
#' network. Where rownames contain unique gene ids and column names are community names
#' @param top_n the number of top genes based on membership scores with higher scores
#' indicating stronger membership in a community
#' @param memb_cut membership threshold for stricter thresholding. Only genes with
#' membership score greater the threshold are used
#' @export
#' @examples
#' \dontrun{
#' memb_list <- list(
#'   GSE39582 = GSE39582_icwgcna$community_membership,
#'   READ = read_icwgcna$community_membership,
#'   COAD = coad_icwgcna$community_membership
#' )
#' ma <- construct_rbh_overlap_based(memb_list)
#' }
construct_rbh_overlap_based <- function(network_membership_list,
                                        top_n = 50,
                                        memb_cut = 0) {
  metaStudies <- names(network_membership_list)
  ns <- sapply(network_membership_list, ncol)
  N <- sum(ns)
  consensus_comms <- unlist(sapply(
    names(network_membership_list),
    function(x) {
      paste0(
        x, "_",
        colnames(network_membership_list[[x]])
      )
    }
  ))
  rbh <- matrix(0, N, N)
  colnames(rbh) <- consensus_comms
  rownames(rbh) <- consensus_comms

  rbh_metrics <- matrix(NA, choose(length(network_membership_list), 2), 5)
  cnt <- 1

  for (i in 1:(length(network_membership_list) - 1))
  {
    for (j in (i + 1):length(network_membership_list))
    {
      rowStart <- sum(ns[(1:i) - 1]) + 1
      rowStop <- rowStart + ns[i] - 1
      colStart <- sum(ns[(1:j) - 1]) + 1
      colStop <- colStart + ns[j] - 1
      tempRBH <- construct_2study_rbh_overlap_based(
        network_membership_list[[i]],
        network_membership_list[[j]],
        top_n = top_n,
        memb_cut = memb_cut
      )

      rbh[rowStart:rowStop, colStart:colStop] <- tempRBH
      message(
        metaStudies[i], ": ",
        ncol(network_membership_list[[i]]),
        " coms, ", metaStudies[j], ": ",
        ncol(network_membership_list[[j]]), " coms, RBH: ",
        sum(apply(tempRBH > 0, 1, sum)), " coms"
      )
      rbh_metrics[cnt, ] <- c(
        metaStudies[i],
        metaStudies[i],
        ncol(network_membership_list[[i]]),
        ncol(network_membership_list[[j]]),
        sum(apply(tempRBH > 0, 1, sum))
      )
      cnt <- cnt + 1
    }
  }

  rbh <- rbh + t(rbh) # contains overlap counts for reciprical best hits
  return(rbh)
}

#' Construct Meta Reciprocal Best Hits based on correlations
#'
#' Iterate through a list of membership matrices and construct pairwise
#' networks for every combination, and output them as a symmetric matrix with
#' labeled dimensions. This can be interpreted as an "adjacency" matrix between
#' all of the different metagenes across all the datasets in the file list. The
#' dimension names will take the original dimension names and concatenate them
#' with the unique file names, such that the dimension names are guaranteed to
#' be unique.
#'
#' @param network_membership_list a list containing community membership scores for each
#' network. Where rownames contain unique gene ids and column names are community names
#' @param lower_quant indicates the quantile for the minimum correlation
#' for the reciprocal best hits we will find.
#' @param upper_quant indicates the quantile for correlations above which
#' ANY metagene pairing will be considered a "hit."
#' @param max_rank represents highest column and row rankings accepted for our
#' reciprocal best hits network. Pure reciprocal best hits uses max_rank of 1.
#' @param abs logical, take absolute values of correlations?
#' @param sparse logical, use a sparse matrix to store network?
#' @param method string, same as [stats::cor()]
#' @param binary logical, indicates whether or not the meta residual best hits
#' matrix should show correlations or simply binary (as.numeric(correlation >
#' 0)) output
#'
#' @return reciprocal best hits matrix between all metagenes across all of the
#' datasets in the specified files, with dimensions named uniquely based on the
#' file column names as well as the file names (to ensure uniqueness).
#' @export
#'
#' @examples
#' \dontrun{
#' memb_list <- list(
#'   GSE39582 = GSE39582_icwgcna$community_membership,
#'   READ = read_icwgcna$community_membership,
#'   COAD = coad_icwgcna$community_membership
#' )
#' ma <- construct_rbh_correlation_based(
#'   memb_list,
#'   upper_quant = .99,
#'   lower_quant = .05,
#'   max_rank = 2
#' )
#' }
construct_rbh_correlation_based <- function(network_membership_list,
                                            lower_quant = 0,
                                            upper_quant = 1.0,
                                            max_rank = 1,
                                            abs = FALSE,
                                            sparse = FALSE,
                                            method = "pearson",
                                            binary = FALSE) {
  metaStudies <- names(network_membership_list)
  ns <- sapply(network_membership_list, ncol)
  N <- sum(ns)
  consensus_comms <- unlist(sapply(
    names(network_membership_list),
    function(x) {
      paste0(
        x, "_",
        colnames(network_membership_list[[x]])
      )
    }
  ))
  rbh <- matrix(0, N, N)
  colnames(rbh) <- consensus_comms
  rownames(rbh) <- consensus_comms

  rbh_metrics <- matrix(NA, choose(length(network_membership_list), 2), 5)
  cnt <- 1

  for (i in 1:(length(network_membership_list) - 1))
  {
    for (j in (i + 1):length(network_membership_list))
    {
      rowStart <- sum(ns[(1:i) - 1]) + 1
      rowStop <- rowStart + ns[i] - 1
      colStart <- sum(ns[(1:j) - 1]) + 1
      colStop <- colStart + ns[j] - 1

      tempRBH <- construct_2study_rbh_correlation_based(
        net_membership_1 = network_membership_list[[i]],
        net_membership_2 = network_membership_list[[j]],
        lower_quant = lower_quant,
        upper_quant = upper_quant,
        max_rank = max_rank,
        abs = abs,
        sparse = sparse,
        method = method
      )

      rbh[rowStart:rowStop, colStart:colStop] <- tempRBH
      message(
        metaStudies[i], ": ",
        ncol(network_membership_list[[i]]),
        " coms, ", metaStudies[j], ": ",
        ncol(network_membership_list[[j]]), " coms, RBH: ",
        sum(apply(tempRBH > 0, 1, sum)), " coms"
      )
      rbh_metrics[cnt, ] <- c(
        metaStudies[i],
        metaStudies[i],
        ncol(network_membership_list[[i]]),
        ncol(network_membership_list[[j]]),
        sum(apply(tempRBH > 0, 1, sum))
      )
      cnt <- cnt + 1
    }
  }
  # contains overlap counts for reciprocal best hits
  rbh <- rbh + t(rbh)
  return(rbh)
}

####
#' RBH Heatmap Creation
#'
#' @param rbh reciprocal best hits matrix (output of
#' [construct_rbh_correlation_based()] or
#' [construct_rbh_overlap_based()])
#' @param memb_list memb_list?
#' @param file_name File name
#' @param width figure width
#' @param height figure height
#'
#' @return pheatmap object
#' @export
#' @examples
#' \dontrun{
#' memb_list <- list(
#'   GSE39582 = GSE39582_icwgcna$community_membership,
#'   READ = read_icwgcna$community_membership,
#'   COAD = coad_icwgcna$community_membership
#' )
#' ma <- construct_rbh_overlap_based(memb_list)
#' plot_rbh(ma, memb_list)
#' }
plot_rbh <- function(rbh, memb_list, file_name = NA, width = 10, height = 8) {
  check_installed("pheatmap")

  ns <- sapply(memb_list, ncol)
  gaps <- sapply(1:length(ns), function(i) {
    sum(ns[1:i])
  })
  anns <- data.frame(Cohorts = gsub("_m.*$", "", rownames(rbh)))
  rownames(anns) <- rownames(rbh)

  pheatmap::pheatmap(
    rbh,
    cluster_rows = FALSE, cluster_cols = FALSE,
    gaps_row = gaps, gaps_col = gaps,
    show_rownames = FALSE, show_colnames = FALSE,
    annotation_row = anns,
    annotation_col = anns,
    color = grDevices::colorRampPalette(c("lightgrey", "blue", "navy"))(50),
    filename = file_name, width = width, height = height
  )
}
