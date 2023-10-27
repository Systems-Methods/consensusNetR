#' Detect Communities in Adjacency/Reciprocal Best Hits Matrix
#'
#' Performs Markov Clustering on the adjacency matrix, giving the
#' cluster/community number that each metagene is associated with, renamed so
#' that the clusters are in order of size (1 being the largest)
#' @param rbh_mat adjacency/reciprocal best hits matrix, as output by
#' [construct_rbh_correlation_based()]
#' @param expansion expansion parameter
#'
#' @return The output of [MCL::mcl()], renamed so that the clusters
#' are in order of size (1 being the largest)
#' @export
#'
#' @examples
#' \dontrun{
#' network_file_dir <- system.file("extdata", package = "consensusNetR")
#' network_file_list <- list.files(network_file_dir, full.names = TRUE)
#' ma <- construct_rbh_correlation_based(
#'   network_file_list = network_file_list,
#'   upper_quant = .99,``
#'   lower_quant = .05, max_rank = 2
#' )
#' comms <- detect_consensus_communities(ma)
#' }
detect_consensus_communities <- function(rbh_mat,
                                        expansion = 2) {

  comms <- MCL::mcl(rbh_mat,
                    addLoops = TRUE,
                    allow1 = FALSE,
                    expansion = expansion)

  if (length(grep('Error', rbh_mat)) > 0) {
    stop('rbh_mat could not be transformed into an equilibrium state matrix. Set expansion parameter to a higher value than 2.')
  }

  cat("Clusters detected\n")
  mapping_list <- as.list(seq_len(length(unique(comms$Cluster))))
  names(mapping_list) <- names(sort(table(comms$Cluster),
                                    decreasing = TRUE))

  comms$Cluster <- as.vector(unlist(mapping_list[as.character(comms$Cluster)]))
  return(comms)
}

