#' Plot individual eigengene distributions
#' @param eigen_list a list of community signatures (eigengenes) for studies,
#' using the output of [calc_consensus_memberships()] as the `membership_matrix`.
#' @param target_study_index target_study_index
#' @param filename File name
#' @param device Device to use. Can either be a device function (e.g. png), or one of "eps", "ps", "tex" (pictex), "pdf", "jpeg", "tiff", "png", "bmp", "svg" or "wmf" (windows only).
#' @param width Plot width in units ("in", "cm", "mm", or "px"). If not supplied, uses the size of current graphics device.
#' @param height Plot height in units ("in", "cm", "mm", or "px"). If not supplied, uses the size of current graphics device.
#' @param dpi Plot resolution. Also accepts a string input: "retina" (320),
#' "print" (300), or "screen" (72). Applies only to raster output types.
#' @return eigengene distributions if `filename` is na, and saved figure is
#' `filename` provided
#' @details
#' In icWGCNA `eigen_list` is created using
#' [icWGCNA::compute_eigengene_matrix()](https://systems-methods.github.io/icWGCNA/reference/compute_eigengene_matrix.html)
#' for each study and then combining all outputs into a list.
#'
#' @export
#' @importFrom rlang .data
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
#' consensus_memb <- get_gene_community_membership(consensus_comms, memb_list, 2)
#'
#' # Need to use icWGCNA for individual eigengenes
#' GSE39582_eigen <- icWGCNA::compute_eigengene_matrix(
#'   ex = GSE39582_df,
#'   membership_matrix = consensus_memb
#' )
#' read_eigen <- icWGCNA::compute_eigengene_matrix(
#'   ex = read_df,
#'   membership_matrix = consensus_memb
#' )
#' coad_eigen <- icWGCNA::compute_eigengene_matrix(
#'   ex = coad_df,
#'   membership_matrix = consensus_memb
#' )
#' eigen_list <- list(GSE39582_eigen, read_eigen, coad_eigen)
#'
#' plot_consensus_eig_dist(eigen_list)
#' }
plot_consensus_eig_dist <- function(eigen_list,
                                    target_study_index = 1,
                                    filename = NA,
                                    device = "png",
                                    width = 12,
                                    height = 10,
                                    dpi = 1000) {
  check_installed(c("ggplot2", "hrbrthemes", "tidyr"))

  temp_eigen <- as.data.frame(t(as.data.frame(lapply(
    normalize_eigengenes(
      eigen_list = eigen_list,
      target_study_index = target_study_index
    ),
    unlist
  ))))
  temp_eigen$Community <- factor(rownames(temp_eigen),
    levels = rownames(temp_eigen)
  )
  temp_eigen <- tidyr::pivot_longer(
    data = temp_eigen, names_to = "sample",
    cols = colnames(temp_eigen)[-ncol(temp_eigen)]
  )

  eigen_dens_plots <- ggplot2::ggplot(
    temp_eigen,
    ggplot2::aes(
      x = .data$value,
      color = .data$Community,
      fill = .data$Community
    )
  ) +
    ggplot2::geom_density(alpha = 0.8) +
    ggplot2::facet_wrap(~Community, scales = "free") +
    hrbrthemes::theme_ipsum() +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::theme(
      legend.position = "none",
      panel.spacing = ggplot2::unit(0.1, "lines"),
      strip.text.x = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    )

  if (is.na(filename)) {
    print(eigen_dens_plots)
  } else {
    ggplot2::ggsave(filename, eigen_dens_plots,
      device = device,
      height = height, width = width, dpi = dpi
    )
  }
}

# quantile normalize of each communities eigen genes based on a target network
normalize_eigengenes <- function(eigen_list, target_study_index = 1) {
  check_installed("aroma.light")

  consensus_comms <- rownames(eigen_list[[target_study_index]])
  qnormed_list <- list()
  for (i in 1:nrow(eigen_list[[target_study_index]])) {
    t_list <- plyr::llply(eigen_list, function(x) {
      return(x[i, ])
    })
    q_normed <- aroma.light::normalizeQuantileRank(t_list, xTarget = sort(t_list[[target_study_index]]))
    qnormed_list[[consensus_comms[i]]] <- q_normed
  }
  return(qnormed_list)
}
