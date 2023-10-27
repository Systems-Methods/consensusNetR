#' Plot individual eigengene distributions
#' @param eigen_list eigen_list
#' @param target_study_index target_study_index
#' @param filename File name
#' @param width figure width
#' @param height figure height
#' @param dpi Plot resolution. Also accepts a string input: "retina" (320),
#' "print" (300), or "screen" (72). Applies only to raster output types.
#' @return eigengene distributions if `filename` is na, and saved figure is
#' `filename` provided
#' @export
#' @importFrom rlang .data
plot_consensus_eig_dist <-  function(eigen_list,
                                     target_study_index = 1,
                                     filename = NA,
                                     width = 12,
                                     height = 10,
                                     dpi = 1000)
{
  needed_packages <- c('ggplot2', 'hrbrthemes')
  missing_packages <- !vapply(needed_packages,
                              FUN = requireNamespace, quietly = TRUE,
                              FUN.VALUE = logical(1))
  if (any(missing_packages)) {
    stop('Must have the following R packages installed for this function: ',
         paste0(names(missing_packages[missing_packages]), collapse = ', '))
  }

  qnormed_list   <- lapply(normalize_eigengenes(eigen_list = eigen_list, target_study_index = target_study_index),unlist)
  temp           <- as.data.frame(qnormed_list); rm(qnormed_list)

  # Plots for eigengene distributions
  temp_eigen           <- as.data.frame(t(temp));
  temp_eigen$Community <- factor(rownames(temp_eigen), levels = rownames(temp_eigen))
  temp_eigen           <- tidyr::pivot_longer(data = temp_eigen,names_to = "sample",
                                              cols = colnames(temp_eigen)[-ncol(temp_eigen)])

  eigen_dens_plots  <- ggplot2::ggplot(
    temp_eigen,
    ggplot2::aes(x = .data$value,
                 color = .data$Community,
                 fill = .data$Community)) +
    ggplot2::geom_density(alpha = 0.8) +
    ggplot2::facet_wrap(~Community,scales = "free") +
    hrbrthemes::theme_ipsum() +
    ggplot2::ylab("") +
    ggplot2::xlab("") +
    ggplot2::theme(
      legend.position = "none",
      panel.spacing = ggplot2::unit(0.1, "lines"),
      strip.text.x = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_blank())

  if (is.na(filename)) {
    eigen_dens_plots
  } else {
    ggplot2::ggsave(filename,eigen_dens_plots,
                    height = height, width = width, dpi = dpi)
  }
}

# quantile normalize of each communities eigen genes based on a target network
normalize_eigengenes <- function(eigen_list, target_study_index=1){
  needed_packages <- c('aroma.light')
  missing_packages <- !vapply(needed_packages,
                              FUN = requireNamespace, quietly = TRUE,
                              FUN.VALUE = logical(1))
  if (any(missing_packages)) {
    stop('Must have the following R packages installed for this function: ',
         paste0(names(missing_packages[missing_packages]), collapse = ', '))
  }

  comms        <- rownames(eigen_list[[target_study_index]])
  qnormed_list <- list()
  for (i in 1:nrow(eigen_list[[target_study_index]])) {
    t_list   <- plyr::llply(eigen_list, function(x){return(x[i,])})
    q_normed <- aroma.light::normalizeQuantileRank(t_list, xTarget = sort(t_list[[target_study_index]]))
    qnormed_list[[comms[i]]] <- q_normed
  }
  return(qnormed_list)
}


