#' Get Consensus Membership of Genes in Metagene Communities
#'
#' @param consensus_comms communities output from [MCL::mcl()],
#' a list containing elements `Cluster` with a vector of the cluster
#' associated with each index.
#' @param network_membership_list a list containing community membership
#' scores for each
#' network. Where rownames contain unique gene ids and column names are community names
#' @param min_studies the minimum threshold of studies/datasets indicating a
#' gene belongs to a community. In a given dataset, a gene "belongs" to a
#' community if the gene's highest loading corresponds to that community.
#' @param include_nonmembers include_nonmembers
#' which do not appear in a community.
#' @param rank_based flag indicating whether to use rank when determing max metagene.
#' Currently only TRUE is supported
#' @param compress indicates whether to drop duplicate community metagenes. Assumes
#' first one is the one to keep (largest or best). Currently only TRUE is supported
#'
#' @return data.frame containing gene ids and their associated most likely
#' cluster(s) in $cluster and $n_studies, representing the number of studies
#' or datasets indicating a gene belongs to the most associated community. In a
#' given dataset, a gene "belongs" to a community if the gene's highest loading
#' corresponds to that community.
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
#' consensus_genes <- get_gene_community_membership(consensus_comms, memb_list, 2)
#' }
get_gene_community_membership <- function(consensus_comms,
                                          network_membership_list,
                                          min_studies = 2,
                                          include_nonmembers = FALSE,
                                          compress = TRUE,
                                          rank_based = TRUE) {
  n_metagenes <- lapply(network_membership_list, ncol)

  gene_cluster_map <- get_gene_map_per_dataset(
    consensus_comms,
    network_membership_list,
    n_metagenes,
    compress = TRUE,
    rank_based = TRUE
  )
  message("Gene cluster map completed")
  get_consensus_mem_from_votes(
    gene_cluster_map,
    min_studies,
    include_nonmembers
  )
}


get_consensus_mem_from_votes <- function(gene_cluster_votes,
                                         min_studies,
                                         include_nonmembers = FALSE) {
  consensus_list <- lapply(
    seq_len(nrow(gene_cluster_votes)),
    function(row_num) {
      gene_id <- row.names(gene_cluster_votes)[row_num]
      counts <- table(gene_cluster_votes[row_num, ])
      clusters <- paste(
        names(counts)[which(counts == max(counts) &
          counts >= min_studies)],
        collapse = ";"
      )
      n_studies <- max(counts)
      if (nchar(clusters) == 0 & !include_nonmembers) {
        return(NULL)
      }
      return(data.frame(
        gene_id = gene_id, cluster = clusters,
        n_studies = n_studies
      ))
    }
  )
  return(do.call(rbind, consensus_list))
}

get_gene_map_per_dataset <- function(consensus_comms,
                                     network_membership_list,
                                     n_metagenes,
                                     compress = TRUE,
                                     rank_based = TRUE) {
  gene_ids <- unique(unlist(lapply(network_membership_list, rownames)))

  ents_clusts <- matrix(NA, length(gene_ids), ncol = length(network_membership_list))
  rownames(ents_clusts) <- gene_ids
  start <- 0
  stop <- 0
  for (i in seq_len(length(network_membership_list))) {
    start <- stop + 1
    stop <- start + n_metagenes[[i]] - 1 # this was missing "- 1" which may have caused many issues
    loc_comm <- consensus_comms$Cluster[start:stop]

    message("loading ", names(network_membership_list)[i])
    network <- network_membership_list[[i]]
    if (compress) {
      # this works under assumption that the first unique metagene is the one to keep.
      keep_inds <- !duplicated(loc_comm)
      loc_comm <- loc_comm[keep_inds]
      network <- network[, keep_inds]
    }
    # For each gene in the i-th network, find the cluster with the highest
    # loading for that gene
    if (rank_based) {
      network <- -apply(-network, 2, rank)
    }
    mod <- apply(network, 1, function(x) {
      which(x == max(x, na.rm = TRUE))[1]
    })
    # For each module in the i-th network, find that module's cluster
    # (these are located in loc_comm) and assign that cluster as the label
    # in the i-th column for that gene in entClusts.
    # This gives us, for each gene, a vote from each dataset for the cluster
    # which that gene is best associated with
    mods <- data.frame(id = row.names(network), mod = loc_comm[mod])
    ents_clusts[, i] <- mods[match(rownames(ents_clusts), mods$id), "mod"]
  }

  ents_clusts <- ents_clusts[apply(!is.na(ents_clusts), 1, sum) > 1, ]

  return(ents_clusts)
}
