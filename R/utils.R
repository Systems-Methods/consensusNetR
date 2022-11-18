# utility function for Fisher's z transform
fisherZ <- function(r){.5*(log(1 + r) - log(1 - r))}
####

#  function for computing cor of cor between study expression data
# takes list of expression datasets in as input
calc_cor_of_cor <- function(ex_list){
  genes <- rownames(ex_list[[1]]); for(i in 2:length(ex_list)){ genes <- genes[genes %in% rownames(ex_list[[i]])] }
  doMC::registerDoMC(parallel::detectCores())
  
  message('Computing intra-study correlations...')
  
  cor_mat <- plyr::llply(ex_list, .fun = function(ex){
    temp  <- as.matrix(ex)
    t_ex  <- temp[intersect(genes,rownames(temp)),]; rm(temp); gc()
    ret   <- as.vector(as.dist(Rfast::cora(t(t_ex)))); rm(t_ex); gc() # Rfast is your friend
    return(ret) 
  }, .parallel = T) 
  rm(ex_list); gc()
  cor_mat <- do.call(cbind, cor_mat)
  
  message('Converting correlations to Z scores...')
  cor_mat     <- DescTools::FisherZ(cor_mat) # transform correlations in to z space in order to compute correlation bewteen datasets.
  
  remove_indices     <- apply(cor_mat, 2, function(i) which(is.na(i)))
  remove_indices_inf <- apply(cor_mat, 2, function(i) which(is.infinite(i)))
  remove_indices     <- unique(c(unlist(remove_indices), unlist(remove_indices_inf)))
  if (length(remove_indices) > 0){ cor_mat <- cor_mat[-remove_indices,] }
  
  cor_of_cors            <- cor(cor_mat, use="pairwise.complete.obs")
  colnames(cor_of_cors)  <- colnames(cor_mat)
  row.names(cor_of_cors) <- colnames(cor_mat)
  rm(cor_mat); gc()
  
  return(cor_of_cors)
}

###

compare_networks <- function(net_memb_1,
                             net_memb_2,
                             K = 75,
                             memb_cut = 0.5,
                             na_flag = "none") {
  K_1        <- min(ncol(net_memb_1), K)
  K_2        <- min(ncol(net_memb_2), K)
  net_memb_1 <- net_memb_1[rownames(net_memb_1) %in% rownames(net_memb_2),1:K_1]
  net_memb_2 <- net_memb_2[match(rownames(net_memb_1), rownames(net_memb_2)),1:K_2]
  rbh        <- compute_2Network_RBH_Overlap_Based(net_memb_1, net_memb_2) 
  overlap    <- sum(rbh > 0)
  
  comms_1    <- apply(net_memb_1, 1,
                      function(x){
                        x[x < memb_cut] <- 0
                        ind <- which(x == max(x))
                        if (length(ind) > 1) {ind <- NA}
                        return(ind)})
  comms_2    <- apply(net_memb_2, 1,
                      function(x){
                        x[x < memb_cut] <- 0;
                        ind <- which(x == max(x))
                        if (length(ind) > 1) {ind <- NA}
                        return(ind)})

  if (na_flag == "none") {
    inds <- rep(TRUE,length(comms_1))
  } else if (na_flag == "both" ) {
    inds <- (!is.na(comms_1)) & !(is.na(comms_2))
  } else if (na_flag == "either") {
    inds <- (!is.na(comms_1)) | (!is.na(comms_2))
  } else {
    stop("na_flag must be \"none\", \"both\", or \"either\"" )
  }

  comms_1    <- comms_1[inds]
  comms_2    <- comms_2[inds]
  net_memb_1 <- net_memb_1[inds,]; net_memb_1[net_memb_1 == 1] <- .99; net_memb_1[net_memb_1 == -1] <- -.99;
  net_memb_2 <- net_memb_2[inds,]; net_memb_2[net_memb_2 == 1] <- .99; net_memb_2[net_memb_2 == -1] <- -.99;
  adjRand    <- mclust::adjustedRandIndex(comms_1,comms_2)

  net_memb_1 <- fisherZ(net_memb_1)
  net_memb_2 <- fisherZ(net_memb_2)
  dist_1     <- stats::as.dist(Rfast::cora(t(net_memb_1)))
  dist_2     <- stats::as.dist(Rfast::cora(t(net_memb_2)))
  cor_mat    <- cbind(fisherZ(unlist(dist_1)), fisherZ(unlist(dist_2)))
  
  remove_indices     <- apply(cor_mat, 2, function(i) which(is.na(i)))
  remove_indices_inf <- apply(cor_mat, 2, function(i) which(is.infinite(i)))
  remove_indices     <- unique(c(unlist(remove_indices), unlist(remove_indices_inf)))
  if (length(remove_indices) > 0){ cor_mat <- cor_mat[-remove_indices,] }
  
  cor_of_cor <- stats::cor(cor_mat[,1], cor_mat[,2],use = "pairwise")

  dend_1     <- stats::hclust(1 - dist_1)
  dend_2     <- stats::hclust(1 - dist_2)
  cor_coph   <- dendextend::cor_cophenetic(dend_1, dend_2)
  
  return(data.frame(adjRand, cor_of_cor, cor_coph, overlap))
}


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
#' @export
#'
compute_2Network_RBH_Overlap_Based <- function(meta_g1,
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




#' Compute RBH for Two Metagene Datasets
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
#' @inheritParams construct_meta_rbh
#'
#' @return Matrix with rows matching the columns in meta_g1 and columns
#' matching the columns of meta_g2, with values 1 indicating two metagenes are a
#' hit, and 0 indicating that they are not.
#' @export
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
#' rbh <- compute_2study_rbh_Correlation_Based(ds1, ds2)
#' }
compute_2study_rbh_Correlation_Based <- function(meta_g1, meta_g2, lower_quant = 0,
                                                 upper_quant = 1.0, max_rank = 1,
                                                 abs = FALSE, sparse = TRUE, method = "pearson") {

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
  return(rbh_network)
}

######## Construct Multi Network Overlap Bases RBH Matrix
#' Construct Meta Reciprocal Best Hits bases on overlaps
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
#' network. Where rownames contain unique gene ids and column names are commuity names
#' @param top_n the number of top genes based on membership scores with higher scores
#' indicating stronger membership in a community
#' @param memb_cut membership threshold for stricter thresholding. Only genes with
#' membership score greater the threshold are used
construct_multi_rbh_overlap_based <- function(network_membership_list,
                                              top_n = 50,
                                              memb_cut = 0) {
  metaStudies <- names(network_membership_list)
  ns          <- sapply(network_membership_list, ncol)
  N           <- sum(ns)
  comms       <- unlist(sapply(names(network_membership_list),
                               function(x){
                                 paste0(x,"_",
                                        colnames(network_membership_list[[x]]))}))
  rbh         <- matrix(0, N,N)
  colnames(rbh) <- comms
  rownames(rbh) <- comms;

  rbh_metrics <- matrix(NA, choose(length(network_membership_list),2), 5);
  cnt         <- 1

  for (i in 1:(length(network_membership_list) - 1))
  {
    for (j in (i + 1):length(network_membership_list))
    {
      rowStart <- sum(ns[(1:i) - 1]) + 1
      rowStop  <- rowStart + ns[i] - 1
      colStart <- sum(ns[(1:j) - 1]) + 1
      colStop  <- colStart + ns[j] - 1
      tempRBH  <- compute_2Network_RBH_Overlap_Based(
        network_membership_list[[i]],
        network_membership_list[[j]],
        top_n = top_n,
        memb_cut = memb_cut
      )

      rbh[rowStart:rowStop, colStart:colStop] <- tempRBH
      message(metaStudies[i],": ",
              ncol(network_membership_list[[i]]),
              " coms, ", metaStudies[j],": ",
              ncol(network_membership_list[[j]]), " coms, RBH: ",
              sum(apply(tempRBH > 0,1,sum))," coms")
      rbh_metrics[cnt, ] <- c(metaStudies[i],
                              metaStudies[i],
                              ncol(network_membership_list[[i]]),
                              ncol(network_membership_list[[j]]),
                              sum(apply(tempRBH > 0,1,sum)))
      cnt <- cnt + 1
    }
  }

  rbh          <- rbh + t(rbh) # contains overlap counts for reciprical best hits
  return(rbh)
}

#### 
plot_rbh <- function(rbh,memb_list, anns, file_name = NA,w=10,h=8){
  ns            <- sapply(memb_list, ncol)
  gaps          <- sapply(1:length(ns),function(i){sum(ns[1:i])})
  #anns          <- data.frame(Cohorts = gsub("_m.*$","",rownames(rbh))); rownames(anns) <- rownames(rbh)
  pheatmap::pheatmap(rbh, cluster_rows = F, cluster_cols = F, 
                     gaps_row = gaps, gaps_col = gaps,
                     show_rownames = F, show_colnames = F, 
                     annotation_row = anns,
                     annotation_col = anns,
                     color=colorRampPalette(c("lightgrey","blue","navy"))(50), filename = file_name , width = w, height = h)
}


#' Construct Meta Reciprocal Best Hits
#'
#' Iterate through a list of files (readable by fread) and construct pairwise
#' networks for every combination, and output them as a symmetric matrix with
#' labeled dimensions. This can be interpreted as an "adjacency" matrix between
#' all of the different metagenes across all the datasets in the file list. The
#' dimension names will take the original dimension names and concatenate them
#' with the unique file names, such that the dimension names are guaranteed to
#' be unique.
#'
#' @param network_file_list a list of paths to files which can be read by
#' data.table::fread. These constitute the metagene datasets. Each file should
#' have some sort of common gene_id as the first column
#' @param network_file_dir a directory in which all files are relevant to the
#' meta reciprocal best hits, all readable by fread. These will be added to the
#' file list. Each file should have some sort of common gene_id as the first
#' column.
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
#' network_file_dir <- system.file("extdata", package = "consensusNetR")
#' network_file_list <- list.files(network_file_dir, full.names = TRUE)
#' ma <- construct_meta_rbh( network_file_list = network_file_list,
#'   upper_quant = .99, lower_quant = .05, max_rank = 2)
#'   }
# Mike/Jimmy, this function needs to become a wrapper function that calls (construct_multi_rbh_overlap_based or construct_multi_rbh_overlap_based based on some flag
# with options "overlap", "pearson" and "spearman" with "overlap" as default. And we need to move away from file lists and just use lists of membership matrices)
construct_meta_rbh <- function(network_file_list = list(),
                               network_file_dir = NULL, lower_quant = 0,
                               upper_quant = 1.0, max_rank = 1,
                               abs = FALSE, sparse = TRUE, method = "pearson",
                               binary = FALSE) {

  ## FileList with full names or directory containing only relevant files
  if ((length(network_file_list) == 0 & is.null(network_file_dir))) {
    stop("construct_meta_rbh requires either a network_file_list or a
         network_file_dir, or both")
  }

  if (!is.null(network_file_dir)) {
    network_file_list <- append(
      network_file_list,
      list.files(network_file_dir, full.names = TRUE)
    )
  }

  ## List of row matrices for the final matrix, representing the rbh of
  ## network_file_list[rownum] with all other datasets in network_file_list,
  ## unless the combination has been analyzed before (The matrix will be
  ## triangular at this stage)
  adj_rows_list <- lapply(
    seq_len(length(network_file_list)),
    function(file1_index) {
      file1_name <- network_file_list[[file1_index]]
      cat(paste("loading", file1_name, "\n"))
      ds1 <- as.matrix(data.table::fread(file1_name))
      row.names(ds1) <- ds1[, 1]
      ds1 <- ds1[, -1]
      colnames(ds1) <- paste0(colnames(ds1),"_",
                              sub(".*/", "", file1_name))
      ds1 <- apply(ds1, c(1,2), as.numeric)

      ## For each dataset, including ds1, onward in the network_file_list, compute
      ## the 2 study network
      row_list <- lapply(
        (file1_index):length(network_file_list),
        function(file2_index) {
          file2_name <- network_file_list[[file2_index]]
          ds2 <- as.matrix(data.table::fread(file2_name))
          row.names(ds2) <- ds2[, 1]
          ds2 <- ds2[, -1]
          colnames(ds2) <- paste0(colnames(ds2),"_",
                                  sub(".*/", "", file2_name))
          ds2 <- apply(ds2, c(1,2), as.numeric)
          return(compute_2study_rbh_Correlation_Based(
            ds1, ds2, lower_quant, upper_quant,
            max_rank, abs, sparse, method
          ))
        }
      )

      ## Nonempty section will be the triangular part of the matrix containing
      ## non-zero values
      if (length(row_list) > 1) {
        nonempty_section <- do.call(cbind, row_list)
      } else {
        nonempty_section <- row_list[[1]]
      }

      if (file1_index == 1) {
        return(nonempty_section)
      }

      ## Create empty section
      empty_sec_list <- lapply(
        1:(file1_index - 1),
        function(file2_index) {
          file2_name <- network_file_list[[file2_index]]
          ds2 <- as.matrix(data.table::fread(file2_name))
          ds2 <- ds2[, -1]
          colnames(ds2) <- paste0(colnames(ds2),"_",
                                  sub(".*/", "", file2_name))
          return(Matrix::Matrix(0,
                                nrow = ncol(ds1),
                                ncol = ncol(ds2),
                                dimnames = list(colnames(ds1), colnames(ds2))
          ))
        })
      empty_section <- do.call(cbind, empty_sec_list)

      return(cbind(empty_section, nonempty_section))
    })

  ## Combine list of row matrices into one large matrix
  rbh_mat <- do.call(rbind, adj_rows_list)
  rbh_mat[lower.tri(rbh_mat)] <- 0
  ## Make matrix symmetrical
  rbh_mat <- rbh_mat + t(as.matrix(rbh_mat * upper.tri(rbh_mat)))
  if (binary) {
    rbh_mat <- Matrix::Matrix(apply(
      rbh_mat, c(1, 2),
      function(x) as.numeric(x > 0)
    ))
  }

  return(rbh_mat)
}

#' Detect Communities in Adjacency/Reciprocal Best Hits Matrix
#'
#' Performs Markov Clustering on the adjacency matrix, giving the
#' cluster/community number that each metagene is associated with, renamed so
#' that the clusters are in order of size (1 being the largest)
#' @param rbh_mat adjacency/reciprocal best hits matrix, as output by
#' [construct_meta_rbh()]
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
#' ma <- construct_meta_rbh(
#'   network_file_list = network_file_list,
#'   upper_quant = .99,``
#'   lower_quant = .05, max_rank = 2
#' )
#' comms <- detect_metagene_communities(ma)
#' }
detect_metagene_communities <- function(rbh_mat,
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


#' Get Consensus Membership of Genes in Metagene Communities
#'
#' @param comms communities output from MCL::mcl, a list containing element
#' $Cluster with a vector of the cluster associated with each index.
#' @param min_studies the minimum threshold of studies/datasets indicating a
#' gene belongs to a community. In a given dataset, a gene "belongs" to a
#' community if the gene's highest loading corresponds to that community.
#' @param include_nonmembers logical indicating whether or not to include genes
#' which do not appear in a community.
#' @param rank_based flag indicating whether to use rank when determing max metagene.
#' Currently only TRUE is supported
#' @param compress indicates whether to drop duplicate community metagenes. Assumes
#' first one is the one to keep (largest or best). Currently only TRUE is supported
#' @inheritParams construct_meta_rbh
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
#' ma <- construct_meta_rbh(
#'   network_file_list = network_file_list,
#'   upper_quant = .99,
#'   lower_quant = .05, max_rank = 2
#' )
#' comms <- detect_metagene_communities(ma)
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
#' ma <- construct_meta_rbh(
#'   network_file_list = network_file_list,
#'   upper_quant = .99,
#'   lower_quant = .05, max_rank = 2
#' )
#' comms <- detect_metagene_communities(ma)
#' compute_mean_meta_genes(comms,ma,network_file_list = network_file_list)
#' }
compute_mean_meta_genes <- function(consensus,
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



#' Identify Top Gene of Communities that are unique (only belong to one community)
#'
#' @param mGenes metat gene matrix rownames being gene IDs and columns being
#' community loadings (i.e. metagenes).
#' @param K number of unique genes to find.
#' @param max_iter maximum nuber of iterations to avoid looping to meaningless top genes
#'
#' @return data.rame with key = community and value = top gene ID. Top genes will only
#' be associated with one module.
#' @export
#'
find_unique_top_genes <- function(mGenes, K=10,
                                  max_iter=10)
{
  converged <- FALSE
  it        <- 1
  while (it <= max_iter & !converged)
  {
    metaGeneRanks   <- t(plyr::aaply(-mGenes,2,rank))
    multiMembInd    <- apply(metaGeneRanks <= K,1,sum) >  1
    nFound <- sum(apply(metaGeneRanks <= K, 1,any) & !multiMembInd)
    message(nFound,"unique genes found with",
            sum(multiMembInd),"non-unique")
    if (sum(multiMembInd) == 0) {converged <- TRUE}
    mGenes <- mGenes[!multiMembInd,]
    it <- it + 1
  }
  metaGeneRanks     <- t(plyr::aaply(-mGenes,2,rank))
  members           <- plyr::ldply(plyr::alply(metaGeneRanks <= K, 2, which),
                                   names)
  rownames(members) <- members[,1]
  members           <- as.data.frame(t(members[,-1]))
  temp              <- tidyr::gather(members)
  ret               <- plyr::ddply(temp, .variables = "value",
                                   .fun = function(x){paste(sort(x$key),
                                                            collapse = ";")})
  ret <- data.frame(key = ret[,2], value = ret[,1])
  ret <- ret[order(ret$key),]
  return(ret)
}

##########
# quantile normalize of each communities eigen genes based on a target network
normalize_eigengenes <- function(eigen_list, target_study_index=1){
  comms        <- rownames(eigen_list[[target_study_index]])
  qnormed_list <- list()
  for(i in 1:nrow(eigen_list[[target_study_index]])){
    t_list   <- plyr::llply(eigen_list, function(x){return(x[i,])})
    q_normed <- aroma.light::normalizeQuantileRank(t_list, xTarget = sort(t_list[[target_study_index]])) 
    qnormed_list[[comms[i]]] <- q_normed
  }
  return(qnormed_list)
}


########
# plot individual eigengene distributions 
plot_consensus_eig_dist <-  function(eigen_list, target_study_index=1, fileName= NULL)
{
  qnormed_list   <- lapply(normalize_eigengenes(eigen_list = eigen_list, target_study_index = target_study_index),unlist)
  temp           <- as.data.frame(qnormed_list); rm(qnormed_list)

  # Plots for eigengene distributions 
  temp_eigen           <- as.data.frame(t(temp));  rm(temp); gc()
  temp_eigen$Community <- factor(rownames(temp_eigen), levels =rownames(temp_eigen))
  temp_eigen           <- tidyr::pivot_longer(data = temp_eigen,names_to = "sample",cols = colnames(temp_eigen)[-ncol(temp_eigen)])
  
  require(ggplot2)
  require(hrbrthemes)
  eigen_dens_plots     <- ggplot(temp_eigen, aes(x=value,color=Community, fill=Community)) + geom_density(alpha=0.8) + facet_wrap(~Community,scales = "free") +  
    hrbrthemes::theme_ipsum() +
    ylab("") + xlab("") + theme(
      legend.position="none", 
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 12),axis.text.x = element_blank(),axis.text.y  = element_blank())
  
  ggsave(fileName,eigen_dens_plots, height = 10, width = 12, dpi = 1000)
}







