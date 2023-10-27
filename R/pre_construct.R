#' Computing cor of cor between study expression data
#'
#' @param ex_list expression datasets
#'
#' @return cor_of_cors
#' @export
#'
calc_cor_of_cor <- function(ex_list){
  needed_packages <- c('doMC','parallel','Rfast')
  missing_packages <- !vapply(needed_packages,
                              FUN = requireNamespace, quietly = TRUE,
                              FUN.VALUE = logical(1))
  if (any(missing_packages)) {
    stop('Must have the following R packages installed for this function: ',
         paste0(names(missing_packages[missing_packages]), collapse = ', '))
  }

  genes <- rownames(ex_list[[1]]);
  for (i in 2:length(ex_list)) {
    genes <- genes[genes %in% rownames(ex_list[[i]])]
  }
  doMC::registerDoMC(parallel::detectCores())

  message(paste('Computing intra-study correlations using', length(genes), "genes"))

  cor_mat <- plyr::llply(ex_list, .fun = function(ex){
    temp  <- as.matrix(ex)
    ret   <- as.vector(stats::as.dist(Rfast::cora(t(temp))));
    return(ret)
  }, .parallel = T)
  cor_mat <- do.call(cbind, cor_mat)

  message('Converting correlations to Z scores...')
  # transform correlations in to z space in order to compute correlation between datasets.
  cor_mat <- DescTools::FisherZ(cor_mat)

  remove_indices     <- apply(cor_mat, 2, function(i) which(is.na(i)))
  remove_indices_inf <- apply(cor_mat, 2, function(i) which(is.infinite(i)))
  remove_indices     <- unique(c(unlist(remove_indices), unlist(remove_indices_inf)))
  if (length(remove_indices) > 0) {
    cor_mat <- cor_mat[-remove_indices,]
  }

  cor_of_cors            <- Rfast::cora(cor_mat)
  colnames(cor_of_cors)  <- colnames(cor_mat)
  row.names(cor_of_cors) <- colnames(cor_mat)
  rm(cor_mat); gc()
  return(cor_of_cors)
}


#' compare_networks
#'
#' @param net_memb_1 net_memb_1
#' @param net_memb_2 net_memb_2
#' @param K K
#' @param memb_cut memb_cut
#' @param na_flag na_flag
#'
#' @return data.frame: adjRand, cor_of_cor, cor_coph, overlap
#' @export
#'
compare_networks <- function(net_memb_1,
                             net_memb_2,
                             K = 75,
                             memb_cut = 0.5,
                             na_flag = c("none", "both", "either")) {
  na_flag <- match.arg(na_flag)
  needed_packages <- c('Rfast')
  missing_packages <- !vapply(needed_packages,
                              FUN = requireNamespace, quietly = TRUE,
                              FUN.VALUE = logical(1))
  if (any(missing_packages)) {
    stop('Must have the following R packages installed for this function: ',
         paste0(names(missing_packages[missing_packages]), collapse = ', '))
  }

  K_1        <- min(ncol(net_memb_1), K)
  K_2        <- min(ncol(net_memb_2), K)
  net_memb_1 <- net_memb_1[rownames(net_memb_1) %in% rownames(net_memb_2),1:K_1]
  net_memb_2 <- net_memb_2[match(rownames(net_memb_1), rownames(net_memb_2)),1:K_2]
  rbh        <- construct_2study_rbh_overlap_based(net_memb_1, net_memb_2)
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

  inds <- switch(
    na_flag,
    none = rep(TRUE, length(comms_1)),
    both = !is.na(comms_1) & !is.na(comms_2),
    either = !is.na(comms_1) | !is.na(comms_2)
  )

  comms_1    <- comms_1[inds]
  comms_2    <- comms_2[inds]
  net_memb_1 <- net_memb_1[inds,]; net_memb_1[net_memb_1 == 1] <- .99;
  net_memb_1[net_memb_1 == -1] <- -.99;
  net_memb_2 <- net_memb_2[inds,]; net_memb_2[net_memb_2 == 1] <- .99;
  net_memb_2[net_memb_2 == -1] <- -.99;
  adjRand    <- mclust::adjustedRandIndex(comms_1,comms_2)

  net_memb_1 <- fisherZ(net_memb_1)
  net_memb_2 <- fisherZ(net_memb_2)
  dist_1     <- stats::as.dist(Rfast::cora(t(net_memb_1)))
  dist_2     <- stats::as.dist(Rfast::cora(t(net_memb_2)))
  cor_mat    <- cbind(fisherZ(unlist(dist_1)), fisherZ(unlist(dist_2)))

  remove_indices     <- apply(cor_mat, 2, function(i) which(is.na(i)))
  remove_indices_inf <- apply(cor_mat, 2, function(i) which(is.infinite(i)))
  remove_indices     <- unique(c(unlist(remove_indices), unlist(remove_indices_inf)))
  if (length(remove_indices) > 0) {
    cor_mat <- cor_mat[-remove_indices,]
  }

  cor_of_cor <- stats::cor(cor_mat[,1], cor_mat[,2],use = "pairwise")

  dend_1     <- stats::hclust(1 - dist_1)
  dend_2     <- stats::hclust(1 - dist_2)
  cor_coph   <- dendextend::cor_cophenetic(dend_1, dend_2)

  return(data.frame(adjRand, cor_of_cor, cor_coph, overlap))
}
