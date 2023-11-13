#' Detect Communities in Adjacency/Reciprocal Best Hits Matrix
#'
#' Performs Markov Clustering on the adjacency matrix, giving the
#' cluster/community number that each metagene is associated with, renamed so
#' that the clusters are in order of size (1 being the largest)
#' @param rbh adjacency/reciprocal best hits matrix, as output by
#' [construct_rbh_correlation_based()]
#' @param expansion expansion parameter
#'
#' @return The output of [MCL::mcl()], renamed so that the clusters
#' are in order of size (1 being the largest)
#' @export
#'
#' @examples
#' \dontrun{
#' # Create list of community_membership object
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
#' consensus_comms <- detect_consensus_communities(ma)
#' }
detect_consensus_communities <- function(rbh,
                                         expansion = 2) {
  consensus_comms <- MCL::mcl(rbh,
    addLoops = TRUE,
    allow1 = FALSE,
    expansion = expansion
  )

  message("Clusters detected")
  mapping_list <- as.list(seq_len(length(unique(consensus_comms$Cluster))))
  names(mapping_list) <- names(sort(table(consensus_comms$Cluster),
    decreasing = TRUE
  ))

  consensus_comms$Cluster <- as.vector(unlist(mapping_list[as.character(consensus_comms$Cluster)]))
  return(consensus_comms)
}
