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
    if (all(is.na(fix[,"ID"]))) stop("SNP IDs are not 'NA' in the VCF file.")
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




