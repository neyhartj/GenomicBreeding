#' Read in SNP marker genotype data
#'
#' @param filename Name of CSV file with marker allele dosage
#' @param geno A data frame of marker allele dosages. If passed, \code{filename} is ignored.
#' @param ploidy The ploidy level e.g. 2, 4, 6, ..., etc.
#' @param phased If TRUE, the input genotype matrix is phased, with alleles coded as positive integers. The name of each
#' individual/clone column is the id and haplotype, concatenated by \code{sep}.
#' @param sep The haplotype concatenator; see \code{phased}.
#' @param inbred Logical. If unphased marker data is passed, should the function assume that all genotypes are inbred?
#'
#' @details
#' The first columns of the file or geno object should be marker, chrom, cM (and optionally bp, if provided). Subsequent columns
#' contain the allele dosage for individuals/clones, coded 0,1,2. While fractional values are allowed for genomewide prediction,
#' these values are rounded when predicting the genetic variance within crosses.
#'
#'
#' @return
#' Variable of class \code{\link{geno}}.
#'
#' @examples
#' # Use data from the PopVar package
#' data("barley_geno", package = "PopVar")
#' data("cranberry_geno", package = "PopVar")
#' data("cranberry_geno_phased", package = "PopVar")
#'
#' # Read in barley marker genotypes assuming no phasing and non-inbred lines.
#' geno_in <- read_geno(geno = barley_geno)
#'
#' # Read in barley marker genotypes assuming completely inbred lines
#' geno_in <- read_geno(geno = barley_geno, inbred = TRUE)
#'
#' # Read in cranberry marker genotype assuming no phasing
#' geno_in <- read_geno(geno = cranberry_geno)
#'
#' # Read in cranberry phased haplotypes from Beagle
#' geno_in <- read_geno(geno = cranberry_geno_phased, phased = TRUE, sep = "_hap")
#'
#' # Simulated polyploidy data
#' geno <- replicate(1000, rbinom(n = 100, size = 4, prob = runif(1)))
#' geno <- t(geno)
#' colnames(geno) <- paste0("id", 1:100)
#' meta <- data.frame(marker = paste0("m", 1:1000), chrom = 1, cM = seq(0, 100, length.out = 1000))
#' geno <- cbind(meta, geno)
#' geno_in <- read_geno(geno = geno, ploidy = 4)
#'
#' @importFrom utils read.csv capture.output
#'
#' @export
#'
read_geno <- function(filename, geno, ploidy = 2L, phased = FALSE, sep = ".", inbred = FALSE) {

  # Error checking
  ploidy <- as.integer(ploidy)
  stopifnot(is.logical(phased))
  stopifnot(is.logical(inbred))
  stopifnot(is.character(sep))

  if (!missing(geno)) {
    data <- geno
  } else {
    if (!endsWith(tolower(filename), "csv")) stop("'filename' is not a CSV.")
    data <- read.csv(file = filename, check.names = FALSE)
  }

  # Extract the map information
  cols_required <- c("marker", "chrom", "cM")
  if (any(!cols_required %in% colnames(data))) {
    stop("The first three columns of 'filename' or 'geno' must be 'marker', 'chrom', and 'cM'.")
  }
  if (colnames(data)[4] == "bp") {
    map <- data[,1:4]
    geno <- as.matrix(data[,-(1:4)])
  } else {
    map <- data[,1:3]
    geno <- as.matrix(data[,-(1:3)])
  }

  m <- nrow(geno)

  if (phased) {
    # No missing data is allowed
    if (any(is.na(geno))) stop("Missing data is not allowed with phased marker genotypes.")

    # Sum haplotypes to get genotypes
    haplo_mat <- geno
    # Unique column names
    id_names <- colnames(haplo_mat)
    id_names_split <- strsplit(x = id_names, split = sep)
    id_names_split <- sapply(X = id_names_split, FUN = "[[", 1)
    id <- unique(id_names_split)

    geno_mat <- matrix(data = as.numeric(NA), nrow = nrow(haplo_mat), ncol = length(id), dimnames = list(row.names(haplo_mat), id))

    # Iterate over unique ids
    for (idx in id) {
      # Find the position of the matching haplotypes; sum the haplotype dosages
      geno_mat[,idx] <- rowSums(haplo_mat[,which(id_names_split == idx)])
    }

    row.names(haplo_mat) <- row.names(geno_mat) <- data[,1]

    # Round the haplotype dosages
    haplo_mat <- round(haplo_mat)

  } else if (inbred) {
    # This functionality only works with diploids
    if (ploidy > 2) stop("Assuming inbred lines works only when ploidy = 2.")

    # Eliminate hets in the genotype matrix and create a haplotype matrix
    geno_mat <- geno
    # Impute with the population mean
    if (any(is.na(geno_mat))) {
      geno_mat1 <- apply(X = geno_mat, MARGIN = 1, FUN = function(snp) {
        snp_mean <- mean(snp, na.rm = TRUE)
        snp[is.na(snp)] <- snp_mean
        snp
      })
      geno_mat <- t(geno_mat1)
    }

    geno_mat <- ifelse(geno_mat >= 1, 2, 0)
    # Expand to a haplotype matrix
    haplo_mat <- matrix(as.numeric(NA), nrow = nrow(geno_mat), ncol = ncol(geno_mat) * 2)

    id <- colnames(geno_mat)
    colnames(haplo_mat) <- paste0(rep(id, each = 2), rep(c("_hapA", "_hapB"), length(id)))
    # Iterate over id
    for (idx in id) {
      haplo_idx <- cbind(geno_mat[,idx], geno_mat[,idx])
      haplo_idx <- haplo_idx / 2
      colnames(haplo_idx) <- paste0(idx, c("_hapA", "_hapB"))
      haplo_mat[,colnames(haplo_idx)] <- haplo_idx
    }

    row.names(haplo_mat) <- row.names(geno_mat) <- data[,1]

  } else {
    haplo_mat <- matrix(data = as.numeric(0), nrow = 0, ncol = 0)
    geno_mat <- geno
    row.names(geno_mat) <- data[,1]

  }

  if (phased) {
    cat("Using phased marker genotypes.\n")
  } else if (inbred) {
    cat("Using inferred phased marker genotypes assuming completely inbred genotypes.\n")
    # Change phased to TRUE if assuming inbred
    phased <- TRUE
  } else {
    cat("Using unphased marker genotypes.\n")
  }

  n <- ncol(geno_mat)
  cat(sub("X",n,"\nNumber of individuals = X\n"))
  m <- nrow(geno_mat)
  cat(sub("X", m,"Number of markers = X\n"))

  geno_mat <- t(geno_mat)
  haplo_mat <- t(haplo_mat)

  output <- new(Class = "geno", ploidy = as.integer(ploidy), map = map, geno.mat = geno_mat, haplo.mat = haplo_mat,
                phased = phased)

  return(output)
}
