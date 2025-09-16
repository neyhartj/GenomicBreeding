#' S4 class of SNP genotype data
#'
#' @slot ploidy The ploidy level.
#' @slot map A data frame of the genetic map information.
#' @slot geno.mat A matrix of marker genotypes.
#' @slot haplo.mat A matrix of haplotypes.
#' @slot phased Whether the marker data was phased.
#' @slot idx The indices of individuals retained after filtering.
#' @slot idy The indices of markers retained after filtering.
#'
geno <- setClass(Class = "geno",
                 slots = c(ploidy = "integer",
                           map = "data.frame",
                           geno.mat = "matrix",
                           haplo.mat = "matrix",
                           phased = "logical",
                           idx = "integer",
                           idy = "integer"))

#' S4 class of microhaplotype genotype data
#'
#' @slot ploidy The ploidy level.
#' @slot map A data frame of the genetic map information.
#' @slot geno.mat A matrix of marker genotypes.
#' @slot haplo.mat A matrix of haplotypes.
#' @slot phased Whether the marker data was phased.
#' @slot idx The indices of individuals retained after filtering.
#' @slot idy The indices of markers retained after filtering.
#'
geno <- setClass(Class = "geno.mh",
                 slots = c(ploidy = "integer",
                           map = "data.frame",
                           geno.mat = "matrix",
                           haplo.mat = "matrix",
                           phased = "logical",
                           idx = "integer",
                           idy = "integer"))

#' S4 class of marker genotype data with linkage disequilibrium information
#'
#' @slot ploidy The ploidy level.
#' @slot map A data frame of the genetic map information.
#' @slot geno.mat A matrix of marker genotypes.
#' @slot haplo.mat A matrix of haplotypes.
#' @slot phased Whether the marker data was phased.
#' @slot idx The indices of individuals retained after filtering.
#' @slot idy The indices of markers retained after filtering.
#' @slot LD A matrix of pairwise marker LD estimates
#'
geno.LD <- setClass(Class = "geno.LD",
                    slots = c(ploidy = "integer",
                              map = "data.frame",
                              geno.mat = "matrix",
                              haplo.mat = "matrix",
                              phased = "logical",
                              idx = "integer",
                              idy = "integer",
                              LD = "Matrix"))


#' S4 class of marker genotype data with genomic relationship matrix information
#'
#' @slot ploidy The ploidy level.
#' @slot map A data frame of the genetic map information.
#' @slot geno.mat A matrix of marker genotypes.
#' @slot haplo.mat A matrix of haplotypes.
#' @slot phased Whether the marker data was phased.
#' @slot idx The indices of individuals retained after filtering.
#' @slot idy The indices of markers retained after filtering.
#' @slot G The additive genomic relationship matrix.
#' @slot D The dominance relationship matrix.
#'
geno.relmat <- setClass(Class = "geno.relmat",
                        slots = c(ploidy = "integer",
                                  map = "data.frame",
                                  geno.mat = "matrix",
                                  haplo.mat = "matrix",
                                  phased = "logical",
                                  idx = "integer",
                                  idy = "integer",
                                  G = "Matrix",
                                  D = "Matrix"))


#' S4 class of a summary of a marker genotype object
#'
#' @slot ploidy The ploidy level.
#' @slot nMar The number of markers.
#' @slot nInd The number of individuals.
#' @slot mac The minor allele count at each locus.
#' @slot mar.missing The proportion of missing data at each marker locus.
#' @slot ind.missing The proportion of missign data of each individual
#'
geno.summary <- setClass(Class = "geno.summary",
                         slots = c(ploidy = "integer",
                                   nMar = "integer",
                                   nInd = "integer",
                                   mac = "numeric",
                                   mar.missing = "numeric",
                                   ind.missing = "numeric"))

