# utility function for Fisher's z transform
fisherZ <- function(r){.5*(log(1 + r) - log(1 - r))}


#' Function for checking needed packages are installed
#'
#' @param packages_needed vector of packages to check
#'
#' @return nothing, but will throw error if packages not installed
check_installed <- function(packages_needed) {
  missing_packages <- !vapply(packages_needed,
                              FUN = is_installed,
                              FUN.VALUE = logical(1))
  if (any(missing_packages)) {
    stop('Must have the following R packages installed for this function: ',
         paste0(names(missing_packages[missing_packages]), collapse = ', '))
  }
}
# Need this wrapper for testing/mocking purposes
is_installed <- function(packages_needed) {
  requireNamespace(packages_needed, quietly = TRUE)
}
