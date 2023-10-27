
<!-- README.md is generated from README.Rmd. Please edit that file -->

# consensusNetR

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

consensusNetR is an R Package for combining networks into a consensus
network based on the work of [Laura Cantini and Andrei
Zinovyev](https://academic.oup.com/bioinformatics/article/35/21/4307/5426054).
In addition to identifying consensus based on correlation of community
meta-genes (loadings or membership scores), we also implement methods
based on overlap.

## Installation

Install the version from BMS BioGit with:

``` r
remotes::install_github(
  repo = "Systems-Methods/consensusNetR"
)
```

or:

``` r
remotes::install_git(
  "https://github.com/Systems-Methods/consensusNetR"
  )
```

# Example Workflow

This example will create a consensus network from three public datasets:
[GSE39582](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39582),
[TCGA COAD](https://portal.gdc.cancer.gov/projects/TCGA-COAD), and [TCGA
READ](https://portal.gdc.cancer.gov/projects/TCGA-READ)

## Downloading data

For GSE39582 we need to convert from Affymetrix Human Genome U133 Plus
2.0 Array to gene symbols, by using the `icWGCNA::gene_mapping()`
function. This matches with the two TCGA datasets already in gene
symbols.

``` r

library(icWGCNA)

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

saveRDS(GSE39582_hugo, 
        file = '/MY_PATH/data/GSE39582_hugo.RDS')


# TCGA READ
UCSCXenaTools::getTCGAdata(
  project = "READ",
  mRNASeq = TRUE, 
  mRNASeqType = "normalized",
  clinical = TRUE, 
  download = TRUE, 
  destdir = "/MY_PATH/data/"
)


# TCGA COAD
UCSCXenaTools::getTCGAdata(
  project = "COAD",
  mRNASeq = TRUE, 
  mRNASeqType = "normalized",
  clinical = TRUE, 
  download = TRUE, 
  destdir = "/MY_PATH/data/"
)
```

## icWGCNA runs

For icWGCNA runs using defaults, except reducing max iterations to 5 for
demonstration purposes. These runs benefit greatly by using multiple
computer cores.

``` r

# GSE39582
GSE39582_hugo <- readRDS(file = '/MY_PATH/data/GSE39582_hugo.RDS')

GSE39582_icwgcna <- icWGCNA::icwgcna(GSE39582_hugo, maxIt = 5)

saveRDS(GSE39582_icwgcna, 
        file = 'MY_PATH/icWGCNA_results/GSE39582_icwgcna.RDS')


# TCGA READ
read_df <- data.table::fread(
  '/MY_PATH/data/TCGA.READ.sampleMap/HiSeqV2.gz', 
  data.table = FALSE) %>% 
  tibble::column_to_rownames('sample')

read_icwgcna <- icWGCNA::icwgcna(read_df, maxIt = 5)
saveRDS(read_icwgcna, 
        file = '/MY_PATH/icWGCNA_runs/read_icwgcna.RDS')


# TCGA COAD
coad_df <- data.table::fread(
  '/MY_PATH/data/TCGA.COAD.sampleMap/HiSeqV2.gz', 
  data.table = FALSE) %>% 
  tibble::column_to_rownames('sample')

coad_icwgcna <- icWGCNA::icwgcna(coad_df, maxIt = 5)
saveRDS(coad_icwgcna, 
        file = '/MY_PATH/icWGCNA_runs/coad_icwgcna.RDS')
```

## Consensus construction

``` r

# Read in data and create list of community_membership object
GSE39582_icwgcna <- readRDS(file = '/MY_PATH/icWGCNA_runs/GSE39582_icwgcna.RDS')
read_icwgcna <- readRDS(file = '/MY_PATH/icWGCNA_runs/read_icwgcna.RDS')
coad_icwgcna <- readRDS(file = '/MY_PATH/icWGCNA_runs/coad_icwgcna.RDS')

memb_list <- list(
  GSE39582 = GSE39582_icwgcna$community_membership,
  READ = read_icwgcna$community_membership,
  COAD = coad_icwgcna$community_membership
)

# Construct Meta Reciprocal Best Hits based on overlaps
rbh <- construct_rbh_overlap_based(memb_list, top_n = 25)

# RBH Heatmap Creation
plot_rbh(rbh = rbh, memb_list = memb_list)
```

<figure>
<img src="man/figures/README-RBH-1.png"
alt="Reciprocal Best Hits Heatmap" />
<figcaption aria-hidden="true">Reciprocal Best Hits Heatmap</figcaption>
</figure>

``` r


# Detect Communities in Adjacency/Reciprocal Best Hits Matrix
comms <- detect_consensus_communities(rbh) 

# Compute the average metagene across studies for each community
consensus_memb  <- calc_consensus_memberships(
  consensus = comms,
  network_membership_list = memb_list
)

# remove the 1st communities which is a miscellaneous comm
consensus_memb  <- consensus_memb[,-1] 

colnames(consensus_memb) <- paste0("mA",colnames(consensus_memb))
write.csv(consensus_memb,
          file = "/MY_PATH/meta_genes.csv", 
          quote = FALSE, 
          row.names = TRUE)
```

## Code of Conduct

Please note that the icWGCNA project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
