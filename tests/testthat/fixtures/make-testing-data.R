# code to make testing data used
library(consensusNetR)
library(tidyverse)
library(icWGCNA)
library(withr)

# GSE39582
GSE39582 <- GEOquery::getGEO("GSE39582")

# creating annotation file for gene mapping to gene symbols
GSE39582_annotation <- GSE39582@featureData@data |>
  dplyr::select(ID, gene_symbol = `Gene Symbol`) |>
  dplyr::mutate(
    gene_symbol = purrr::map(
      gene_symbol, ~ stringr::str_split(.x, " /// ")[[1]]
    )
  ) %>%
  tidyr::unnest(gene_symbol)

GSE39582_hugo <- icWGCNA::gene_mapping(
  GSE39582@assayData$exprs,
  GSE39582_annotation,
  compress_fun = "highest_mean",
  compress_trans = "log_exp"
)
# GSE39582_hugo=readRDS('~/Downloads/GSE39582_hugo.RDS')

# TCGA READ
UCSCXenaTools::getTCGAdata(
  project = "READ",
  mRNASeq = TRUE,
  mRNASeqType = "normalized",
  clinical = TRUE,
  download = TRUE,
  destdir = "/MY_PATH/data/"
)
read_df <- data.table::fread(
  '~/Downloads/TCGA.READ.sampleMap/HiSeqV2.gz',
  data.table = FALSE) %>%
  tibble::column_to_rownames('sample')


# TCGA COAD
UCSCXenaTools::getTCGAdata(
  project = "COAD",
  mRNASeq = TRUE,
  mRNASeqType = "normalized",
  clinical = TRUE,
  download = TRUE,
  destdir = "/MY_PATH/data/"
)

coad_df <- data.table::fread(
  '~/Downloads/TCGA.COAD.sampleMap/HiSeqV2.gz',
  data.table = FALSE) %>%
  tibble::column_to_rownames('sample')


# Read in data and create list of community_membership object
GSE39582_icwgcna <- icWGCNA::icwgcna(GSE39582_hugo, maxIt = 5)
read_icwgcna <- icWGCNA::icwgcna(read_df, maxIt = 5)
coad_icwgcna <- icWGCNA::icwgcna(coad_df, maxIt = 5)

# GSE39582_icwgcna <- readRDS(file = '~/Downloads/GSE39582_icwgcna.RDS')
# read_icwgcna <- readRDS(file = '~/Downloads/read_icwgcna.RDS')
# coad_icwgcna <- readRDS(file = '~/Downloads/coad_icwgcna.RDS')

# Need to dramatically subset for testing
all_genes <- unique(c(
  rownames(GSE39582_icwgcna$community_membership),
  rownames(read_icwgcna$community_membership),
  rownames(coad_icwgcna$community_membership)
  ))
withr::with_seed(
  523757542,
  {
    gene_subset <- sort(sample(all_genes, 500))
    GSE39582_mem <- GSE39582_icwgcna$community_membership[gene_subset,] |>
      na.omit()
    read_mem <- read_icwgcna$community_membership[gene_subset,] |>
      na.omit()
    coad_mem <- coad_icwgcna$community_membership[gene_subset,] |>
      na.omit()
  }
)

memb_list <- list(
  GSE39582 = GSE39582_mem,
  READ = read_mem,
  COAD = coad_mem
)

saveRDS(memb_list,
        file = testthat::test_path("fixtures","memb_list.rds"))

rbh <- suppressMessages(construct_rbh_overlap_based(memb_list))
consensus_comms <- suppressMessages(detect_consensus_communities(rbh))

saveRDS(consensus_comms,
        file = testthat::test_path("fixtures","consensus_comms.rds"))

# Getting ex subset

GSE39582_sub <- GSE39582_hugo[gene_subset, ] |>
  na.omit()
read_sub <- read_df[gene_subset, ] |>
  na.omit()
coad_sub <- coad_df[gene_subset, ] |>
  na.omit()

ex_list <- list(
  GSE39582 = GSE39582_sub,
  READ = read_sub,
  COAD = coad_sub
)

saveRDS(ex_list,
        file = testthat::test_path("fixtures","ex_list.rds"))



eigen_list <- list(
  GSE39582 = GSE39582_icwgcna$community_signature,
  READ = read_icwgcna$community_signature,
  COAD = coad_icwgcna$community_signature
)

saveRDS(eigen_list,
        file = testthat::test_path("fixtures","eigen_list.rds"))



