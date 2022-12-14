% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{construct_meta_rbh}
\alias{construct_meta_rbh}
\title{Construct Meta Reciprocal Best Hits}
\usage{
construct_meta_rbh(
  network_file_list = list(),
  network_file_dir = NULL,
  lower_quant = 0,
  upper_quant = 1,
  max_rank = 1,
  abs = FALSE,
  sparse = TRUE,
  method = "pearson",
  binary = FALSE
)
}
\arguments{
\item{network_file_list}{a list of paths to files which can be read by
data.table::fread. These constitute the metagene datasets. Each file should
have some sort of common gene_id as the first column}

\item{network_file_dir}{a directory in which all files are relevant to the
meta reciprocal best hits, all readable by fread. These will be added to the
file list. Each file should have some sort of common gene_id as the first
column.}

\item{lower_quant}{indicates the quantile for the minimum correlation
for the reciprocal best hits we will find.}

\item{upper_quant}{indicates the quantile for correlations above which
ANY metagene pairing will be considered a "hit."}

\item{max_rank}{represents highest column and row rankings accepted for our
reciprocal best hits network. Pure reciprocal best hits uses max_rank of 1.}

\item{abs}{logical, take absolute values of correlations?}

\item{sparse}{logical, use a sparse matrix to store network?}

\item{method}{string, same as \code{\link[stats:cor]{stats::cor()}}}

\item{binary}{logical, indicates whether or not the meta residual best hits
matrix should show correlations or simply binary (as.numeric(correlation >
0)) output}
}
\value{
reciprocal best hits matrix between all metagenes across all of the
datasets in the specified files, with dimensions named uniquely based on the
file column names as well as the file names (to ensure uniqueness).
}
\description{
Iterate through a list of files (readable by fread) and construct pairwise
networks for every combination, and output them as a symmetric matrix with
labeled dimensions. This can be interpreted as an "adjacency" matrix between
all of the different metagenes across all the datasets in the file list. The
dimension names will take the original dimension names and concatenate them
with the unique file names, such that the dimension names are guaranteed to
be unique.
}
\examples{
\dontrun{
network_file_dir <- system.file("extdata", package = "consensusNetR")
network_file_list <- list.files(network_file_dir, full.names = TRUE)
ma <- construct_meta_rbh( network_file_list = network_file_list,
  upper_quant = .99, lower_quant = .05, max_rank = 2)
  }
}
