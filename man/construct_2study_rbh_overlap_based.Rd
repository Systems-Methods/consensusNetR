% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rbh_2study.R
\name{construct_2study_rbh_overlap_based}
\alias{construct_2study_rbh_overlap_based}
\title{Compute RBH (reciprocal best hits) for Two Networks based on Overlap of top genes}
\usage{
construct_2study_rbh_overlap_based(
  net_membership_1,
  net_membership_2,
  top_n = 50,
  memb_cut = 0
)
}
\arguments{
\item{net_membership_1}{community membership matrix from first network.}

\item{net_membership_2}{community membership matrix from second network.}

\item{top_n}{the number of top genes based on membership scores with higher scores
indicating stronger membership in a community}

\item{memb_cut}{membership threshold for stricter thresholding. Only genes with
membership score greater the threshold are used}
}
\value{
Matrix with rows matching the columns in net_membership_1 and columns
matching the columns of net_membership_2, with values 1 indicating two metagenes are a
hit, and 0 indicating that they are not.
}
\description{
Compute RBH (reciprocal best hits) for Two Networks based on Overlap of top
community genes. They must have row names containing some sort of common geneID
vocabulary, so that they can be intersected to get common rows.
}
\examples{
\dontrun{
rbh <- construct_2study_rbh_overlap_based(
  GSE39582_icwgcna$community_membership,
  read_icwgcna$community_membership
)
}
}
\keyword{internal}
