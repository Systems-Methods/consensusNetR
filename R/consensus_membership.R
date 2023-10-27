#' Compute the average metagene across studies for each community
#'
#' @param consensus communities output from [MCL::mcl()], a list containing element
#' $Cluster with a vector of the cluster associated with each index.
#' @param network_membership_list a list containing community membership scores for each
#' network. Where rownames contain unique gene ids and column names are community names
#' @param weights weights for weighted average of based on study attributes,
#' not yet implemented
#' @param gene_cohort_N gene_cohort_N
#' @param compressIntra indicates how to deal with multiple metagenes belonging
#' to the same community within one study-level network
#'
#' @return matrix with average loadings (metagenes) across studies included from
#' network_file_list
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
#' calc_consensus_memberships(comms,ma,network_file_list = network_file_list)
#' }
calc_consensus_memberships <- function(consensus,
                                 network_membership_list,
                                 weights = NULL,
                                 gene_cohort_N = 3,
                                 compressIntra = "first")
{

  n_metagenes <- lapply(network_membership_list, ncol)
  u_genes     <-  table(unlist(lapply(network_membership_list, rownames)))
  u_genes     <- sort(names(u_genes[u_genes >= gene_cohort_N]))

  # collect gene loadings to compute average gene loadings (actually kme's) across each dataset
  clusterLoadings   <- list()

  for (c in as.character(sort(unique(consensus$Cluster)))) {
    clusterLoadings[[c]] <- matrix(NA,length(u_genes),
                                   length(network_membership_list))
  }

  start <- 0
  stop  <- 0
  for (i in seq_len(length(network_membership_list)))
  {
    start    <- stop + 1
    stop     <- start + n_metagenes[[i]] - 1
    loc_comm <- consensus$Cluster[start:stop]

    cat(paste("loading", names(network_membership_list)[i], "\n"))
    network <- as.matrix(network_membership_list[[i]])
    mods    <- consensus$Cluster[start:stop]

    for (c in as.character(sort(unique(mods))))
    {
      mod_cnt <- sum(mods == c)
      if (mod_cnt == 0) {
        x <-  matrix(NA,nrow(network),1)
      }
      else if (mod_cnt > 1)
      {
        tMi <- network[,mods == c];
        if (compressIntra == "first") {
          x <- tMi[,1]
        }
        else{x <- apply(tMi,1, mean);}
      }else
      {
        x <- network[,mods == c];
      }

      clusterLoadings[[as.character(c)]][,i] <- x[match(u_genes,
                                                        rownames(network))]
    }
  }

  meanGeneLoadings <- t(plyr::laply(clusterLoadings,
                                    .fun = compressMetaGenes,
                                    method = "none", w = weights))
  rownames(meanGeneLoadings) <- u_genes
  colnames(meanGeneLoadings) <- names(clusterLoadings)
  return(meanGeneLoadings)
}


compressMetaGenes <- function(y,method = "none",
                              w=NULL,
                              dynamic_thresh = .65){
  if (is.null(w)) {
    ret <- apply(y,1,mean,na.rm = TRUE);
  } else if (length(w) == 1 & w == "Dynamic") {
    weights <- apply(y,2,function(x){
      x <- sort(x,decreasing = TRUE)
      return(mean(x[1:5]))
    })
    message(weights)

    weights <- weights - dynamic_thresh
    weights[weights < 0 | is.na(weights)] <- 0

    message(weights)
    ret     <- apply(y, 1, stats::weighted.mean,
                     w = weights, na.rm = TRUE)
  }else {
    ret <- apply(y, 1, stats::weighted.mean,
                 w = w, na.rm = TRUE)
  }

  ret[is.na(ret)] <- mean(ret,na.rm = TRUE)
  return(ret)
}
