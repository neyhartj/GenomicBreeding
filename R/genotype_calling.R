#' Convert a VCF to read count matrices for updog
#'
#' @param vcf.in A \code{vcfR} object.
#' @param variants.keep A vector of variant IDs to keep for output. If \code{NULL}, all variants are kept.
#' @param keep.multiallelic If \code{TRUE}, multi-allelic SNPs are kept, but only the most frequent alternate allele is retained. If \code{FALSE}, multi-allelic SNPs are removed.
#'
#' @return A list with two matrices: \code{refmat} (the reference allele read count matrix) and \code{sizemat} (the total allele read count matrix).
#'
#' @import vcfR
#'
#' @export
#'
vcf2updog <- function(vcf.in, variants.keep = NULL, keep.multiallelic = FALSE) {

  stopifnot(inherits(vcf.in, "vcfR"))
  if (!is.null(variants.keep)) {
    stopifnot(inherits(variants.keep, c("numeric", "character")))
  }
  stopifnot(is.logical(keep.multiallelic))

  # Get fix
  fix <- getFIX(x = vcf.in)

  # Subset variants
  if (is.null(variants.keep)) {
    idx <- seq_len(nrow(vcf.in))
  } else if (is.numeric(variants.keep)) {
    idx <- variants.keep
  } else {
    if (all(is.na(fix[,"ID"]))) stop("SNP IDs are all 'NA' in the VCF file.")
    idx <- which(fix[,"ID"] %in% variants.keep)
  }

  # Subset the vcf
  vcf.in <- vcf.in[idx, ]

  # Now deal with multi-allelic sites
  if (!keep.multiallelic) {
    # Just get bi-alleleic sites
    idx <- which(is.biallelic(x = vcf.in))
    vcf.in <- vcf.in[idx, ]

    # Get the allele depth matrix
    ad <- extract.gt(x = vcf.in, element = "AD", IDtoRowNames = TRUE)
    # Get the reference depth
    refmat <- masplit(myMat = ad, delim = ",", count = 0, record = 1, sort = 0)
    # Get the alternate depth
    altmat <- masplit(myMat = ad, delim = ",", count = 0, record = 2, sort = 0)
    # Compute the size mat
    sizemat <- refmat + altmat

  } else {
    # For multiallelic sites, use the reference allele and the alternate allele at greatest frequency
    idx_bi <- which(is.biallelic(x = vcf.in))
    idx_multi <- setdiff(seq_len(nrow(vcf.in)), idx_bi)

    # Get ad element
    ad <- extract.gt(x = vcf.in, element = "AD", IDtoRowNames = TRUE)
    # Get reference depth
    refmat <- masplit(myMat = ad, delim = ",", count = 0, record = 1, sort = 0)
    # Build the altmat
    altmat <- masplit(myMat = ad, delim = ",", count = 0, record = 2, sort = 0)
    # Set multi-allelic to 0
    altmat[idx_multi, ] <- 0

    # Iterate over multi-allelic SNPs
    for (i in idx_multi) {
      # Get alt matrix for alleles 2, 3, 4
      alt_multi_list_i <- sapply(X = 2:4, FUN = function(j) masplit(myMat = ad[i,, drop = FALSE], delim = ",", count = 0, record = j, sort = 0))
      # Set NAs to 0
      alt_multi_list_i[is.na(alt_multi_list_i)] <- 0
      # Which allele has the highest number of >0 reads
      which_alt_allele <- which.max(colSums(alt_multi_list_i > 0))
      # Extract the alt matrix for this allele
      altmat_i <- masplit(myMat = ad[i,, drop = FALSE], delim = ",", count = 0, record = which_alt_allele + 1, sort = 0)
      # add this to the altmat
      altmat_i[i, ] <- altmat_i[1,]

    } # close the loop

    # Compute the size mat
    sizemat <- refmat + altmat

  } # Close the logical statement

  # Export the matrices
  out <- list(refmat = refmat, sizemat = sizemat)

  return(out)

}


#' Convert a MADC file to objects for genotyping calling in updog
#'
#' @param madc A data.frame of MADC data (read from a .csv).
#' @param botloci A character vector indicating the names of loci for which the sequence is on the bottom strand with respect to the reference genome.
#' See \code{\link[polyRAD]{readDArTag}}.
#' @param target.only Logical. Should only target reference and alternate tags be retained? If TRUE, the output data is bi-allelic SNPs only; if FALSE,
#' the output data is multi-allelic loci.
#'
#' @details
#'
#' This function uses polyRAD to import the microhaplotype data, but converts it
#' to updog format for genotyping calling
#'
#' @importFrom polyRAD readDArTag
#'
#' @export
#'
#'
madc2updog <- function(madc, botloci, target.only = FALSE) {

  if (!missing(madc)) {
    stopifnot(is.data.frame(madc))
    madc <- as.data.frame(madc)
  } else if (!missing(madc.file)) {
    stopifnot(file.exists(madc.file))
    madc <- read.csv(file = madc.file, as.is = TRUE, check.names = FALSE)
  } else {
    stop("One of 'madc' or 'madc.file' must be passed.")
  }

  # Verify the names of the first three columns
  stopifnot(all(names(madc)[1:3] == c("AlleleID", "CloneID", "AlleleSequence")))
  stopifnot(is.logical(target.only))
  stopifnot(is.character(botloci))

  # Summarize and print
  loci_names <- unique(madc$CloneID)
  n_loci <- length(loci_names)
  n_alleles <- nrow(madc)
  n_samples <- ncol(madc) - 3
  n_alleles_per_loci <- sapply(split(madc$AlleleID, madc$CloneID), length)

  # Write the madc file as a temporary file
  tmp_dir <- "tmp/"
  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir, showWarnings = FALSE)
  }
  tmp_file <- file.path(tmp_dir, "tmp_madc.csv")
  write.csv(x = madc, file = tmp_file, row.names = FALSE)

  # If target only, extract the reference and alternate allele depths
  if (target.only) {
    # Index of the target reference and alternate counts
    ref_alt_idx <- grep(pattern = "Ref_[0]{1,}1|Alt_[0]{1,}2", x = madc$AlleleID)
    # Subset
    madc_target <- madc[ref_alt_idx, , drop = FALSE]

    # Create a matrix for refmat and sizemat
    refmat <- sizemat <- matrix(as.numeric(NA), nrow = n_loci, ncol = n_samples, dimnames = list(loci_names, names(madc)[-1:-3]))
    # Fill in the matrix
    ref_df <- madc_target[grep(pattern = "Ref", x = madc_target$AlleleID), ]
    row.names(ref_df) <- ref_df$CloneID
    refmat <- as.matrix(ref_df[-1:-3])

    alt_df <- madc_target[grep(pattern = "Alt", x = madc_target$AlleleID), ]
    row.names(alt_df) <- alt_df$CloneID
    altmat <- as.matrix(alt_df[-1:-3])

    if (nrow(refmat) != nrow(altmat)) stop("The size of the reference matrix is not the same as the alternate matrix. Please check the MADC file.")

    sizemat <- refmat + altmat

    # Export the matrices
    out <- list(refmat = refmat, sizemat = sizemat, allele2loc = seq_len(nrow(sizemat)))
    # Add class
    class(out) <- c(class(out), "updog.biallelic")


  } else {
    # Read into polyRAD
    dat <- readDArTag(file = tmp_file, botloci = botloci, trim.sample.names = "")

    # Get the allele depth matrix
    allele_depth_mat <- t(dat$alleleDepth)
    # Get the anti-allele depth matrix
    not_allele_depth_mat <- t(dat$antiAlleleDepth)
    # The size mat the sum of allele depth and anti-allele depth
    allele_size_mat <- allele_depth_mat + not_allele_depth_mat

    # Replace sample names
    colnames(allele_depth_mat) <- colnames(allele_size_mat) <- names(madc)[-1:-3]

    # Export the matrices
    out <- list(refmat = allele_depth_mat, sizemat = allele_size_mat, allele2loc = dat$alleles2loc)
    # Add class
    class(out) <- c(class(out), "updog.multiallelic")

  }

  # Delete the temporary file and directory
  unlink(x = tmp_file)
  unlink(x = tmp_dir, force = TRUE, recursive = TRUE)
  return(out)

}


#' Call marker genotypes using updog
#'
#' @description
#' This is a simple wrapper function around functions from \code{\link[updog]{updog}} but
#' is made aware of multi-allelic markers that are not handled natively by updog.
#'
#' @param x A list of two matrices from the function \code{\link{madc2updog}}.
#' @param ploidy The ploidy of all individuals in the population.
#' @param model See \code{\link[updog]{multidog}}.
#' @param n.cores Number of CPU cores for parallel processing.
#'
#' @export
#'
run_updog <- function(x, ploidy = 2, model = "norm", n.cores = 1) {
  # Verify class
  stopifnot(inherits(x, c("updog.multiallelic", "updog.biallelic")))

  # Run the matrices through updog

  # ## TESTING ##
  # idx <- which(x$allele2loc %in% c(1:10))
  # updog_out <- multidog(refmat = x$refmat[idx, ], sizemat = x$sizemat[idx, ], ploidy = ploidy, model = model,
  #                       nc = n.cores)
  # ## TESTING ##

  updog_out <- multidog(refmat = x$refmat, sizemat = x$sizemat, ploidy = ploidy, model = model,
                        nc = n.cores)


  # Parse the results if multi-allelic
  updog_inddf <- updog_out$inddf


}



























