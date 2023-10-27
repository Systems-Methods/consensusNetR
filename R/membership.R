#' Get Consensus Membership of Genes in Metagene Communities
#'
#' @param comms communities output from MCL::mcl, a list containing element
#' $Cluster with a vector of the cluster associated with each index.
#' @param min_studies the minimum threshold of studies/datasets indicating a
#' gene belongs to a community. In a given dataset, a gene "belongs" to a
#' community if the gene's highest loading corresponds to that community.
#' @param network_file_list network_file_list
#' @param network_file_dir network_file_dir
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
#' network_file_dir <- system.file("extdata", package = "consensusNetR")
#' network_file_list <- list.files(network_file_dir, full.names = TRUE)
#' ma <- construct_rbh_correlation_based(
#'   network_file_list = network_file_list,
#'   upper_quant = .99,
#'   lower_quant = .05, max_rank = 2
#' )
#' comms <- detect_consensus_communities(ma)
#' consensus <- get_gene_community_membership(comms, 2,
#'   network_file_dir = network_file_dir
#' )
#' }
get_gene_community_membership <- function(comms, min_studies,
                                          network_file_list = list(),
                                          network_file_dir = NULL,
                                          include_nonmembers = FALSE,
                                          compress = TRUE,
                                          rank_based = TRUE) {
  if ((length(network_file_list) == 0 & is.null(network_file_dir))) {
    stop("get_gene_community_membership requires either a network_file_list or a
         network_file_dir, or both")
  }

  if (!is.null(network_file_dir)) {
    network_file_list <- append(
      network_file_list,
      list.files(network_file_dir, full.names = TRUE)
    )
  }

  n_metagenes <- lapply(
    network_file_list,
    function(x) {
      ncol(as.matrix(data.table::fread(x))) - 1
    }
  )
  gene_cluster_map <- get_gene_map_per_dataset(
    comms, n_metagenes,
    network_file_list,
    compress = TRUE,
    rank_based = TRUE
  )
  cat("Gene cluster map completed\n")
  return(get_consensus_mem_from_votes(gene_cluster_map, min_studies,
                                      include_nonmembers))
}


get_consensus_mem_from_votes <- function(gene_cluster_votes, min_studies,
                                         include_nonmembers = FALSE) {
  consensus_list <- lapply(seq_len(nrow(gene_cluster_votes)), function(row_num) {
    gene_id <- row.names(gene_cluster_votes)[row_num]
    counts <- table(gene_cluster_votes[row_num, ])
    clusters <- paste(names(counts)[which(counts == max(counts) &
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
  })
  return(do.call(rbind, consensus_list))
}

get_gene_map_per_dataset <- function(community, n_metagenes,
                                     network_file_list = list(),
                                     network_file_dir = NULL, compress = TRUE,
                                     rank_based = TRUE) {
  ## Takes the output of MCL::mcl
  if ((length(network_file_list) == 0 & is.null(network_file_dir))) {
    stop("get_gene_map_per_dataset requires either a network_file_list or a
         network_file_dir, or both")
  }

  if (!is.null(network_file_dir)) {
    network_file_list <- append(
      network_file_list,
      list.files(network_file_dir, full.names = TRUE)
    )
  }

  if (length(n_metagenes) == 1) {
    n_metagenes <- replicate(length(network_file_list), n_metagenes)
  } else if (length(n_metagenes) != length(network_file_list)) {
    stop("length of n_metagenes must either be 1 or equal to the number
         of files in networkFolder + network_file_list")
  }
  gene_ids <- character()
  for (i in seq_len(length(network_file_list))) {
    new_genes <- unlist(data.table::fread(network_file_list[[1]])[, 1])
    gene_ids <- c(gene_ids, new_genes[sapply(
      new_genes,
      function(x) {
        return(!(x %in% gene_ids))
      }
    )])
  }

  ents_clusts <- matrix(NA, length(gene_ids), ncol = length(network_file_list))
  rownames(ents_clusts) <- gene_ids
  start <- 0
  stop <- 0
  for (i in seq_len(length(network_file_list))) {
    start <- stop + 1
    stop <- start + n_metagenes[[i]] - 1 # this was missing "- 1" which may have caused many issues
    loc_comm <- community$Cluster[start:stop]

    cat(paste("loading", network_file_list[[i]], "\n"))
    network <- as.matrix(data.table::fread(network_file_list[[i]]))
    row.names(network) <- network[, 1]
    network <- network[, -1]
    if (compress) #this works under assumption that the first unique metagene is the one to keep.
    {
      keep_inds <- !duplicated(loc_comm)
      loc_comm  <- loc_comm[keep_inds]
      network   <- network[,keep_inds]
    }
    # For each gene in the i-th network, find the cluster with the highest
    # loading for that gene
    if (rank_based) {network <- -apply(-network,2,rank)}
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
