#' Compute the average metagene across studies for each community
#'
#' @param consensus_comms communities output from [MCL::mcl()],
#' a list containing elements `Cluster` with a vector of the cluster
#' associated with each index.
#' @param network_membership_list a list containing community membership scores for each
#' network. Where rownames contain unique gene ids and column names are community names
#' @param weights weights for weighted average of based on study attributes
#' @param gene_cohort_N gene_cohort_N
#' @param compressIntra indicates how to deal with multiple metagenes belonging
#' to the same community within one study-level network
#' @param remove_misc_comm should the miscellaneous community be removed
#' (1st community)
#' @param comm_prefix the prefix to add to community names
#'
#' @return matrix with average loadings (metagenes) across studies included from
#' network_membership_list
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
#' consensus_comms <- detect_consensus_communities(ma)
#' consensus_memb <- calc_consensus_memberships(consensus_comms, memb_list)
#' }
calc_consensus_memberships <- function(consensus_comms,
                                       network_membership_list,
                                       weights = NULL,
                                       gene_cohort_N = 3,
                                       compressIntra = c("first", "mean"),
                                       remove_misc_comm = TRUE,
                                       comm_prefix = "mA") {
  compressIntra <- match.arg(compressIntra)
  n_metagenes <- lapply(network_membership_list, ncol)
  u_genes <- table(unlist(lapply(network_membership_list, rownames)))
  u_genes <- sort(names(u_genes[u_genes >= gene_cohort_N]))

  # collect gene loadings to compute average gene loadings (actually kme's) across each dataset
  clusterLoadings <- list()

  for (c in as.character(sort(unique(consensus_comms$Cluster)))) {
    clusterLoadings[[c]] <- matrix(
      NA, length(u_genes),
      length(network_membership_list)
    )
  }

  start <- 0
  stop <- 0
  for (i in seq_len(length(network_membership_list))) {
    start <- stop + 1
    stop <- start + n_metagenes[[i]] - 1

    message("loading ", names(network_membership_list)[i])
    network <- as.matrix(network_membership_list[[i]])
    mods <- consensus_comms$Cluster[start:stop]

    for (c in as.character(sort(unique(mods)))) {
      tMi <- network[, mods == c, drop = FALSE]
      if (compressIntra == "first") {
        x <- tMi[, 1]
      } else {
        x <- apply(tMi, 1, mean)
      }
      clusterLoadings[[as.character(c)]][, i] <-
        x[match(u_genes, rownames(network))]
    }
  }

  meanGeneLoadings <- do.call(
    cbind,
    lapply(
      clusterLoadings,
      compressMetaGenes,
      method = "none",
      w = weights
    )
  )

  rownames(meanGeneLoadings) <- u_genes
  colnames(meanGeneLoadings) <- paste0(comm_prefix, names(clusterLoadings))
  if (remove_misc_comm) {
    meanGeneLoadings <- meanGeneLoadings[, -1]
  }
  return(meanGeneLoadings)
}


compressMetaGenes <- function(y, method = "none",
                              w = NULL,
                              dynamic_thresh = .65) {
  if (is.null(w)) {
    ret <- apply(y, 1, mean, na.rm = TRUE)
  } else if (length(w) == 1 && w == "Dynamic") {
    weights <- apply(y, 2, function(x) {
      x <- sort(x, decreasing = TRUE)
      return(mean(x[1:5]))
    })
    message(weights)

    weights <- weights - dynamic_thresh
    weights[weights < 0 | is.na(weights)] <- 0

    message(weights)
    ret <- apply(y, 1, stats::weighted.mean,
      w = weights, na.rm = TRUE
    )
  } else {
    ret <- apply(y, 1, stats::weighted.mean,
      w = w, na.rm = TRUE
    )
  }

  ret[is.na(ret)] <- mean(ret, na.rm = TRUE)
  return(ret)
}
