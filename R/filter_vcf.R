#' Filter genotype calls in a VCF object
#'
#' @description
#' Applies filters to genotype calls and sites in a VCF file. The function is designed to apply
#' allele depth filters to heterozygous genotype calls.
#'
#'
#' @param x A \code{vcfR} object.
#' @param minAD The minimum allele depth to retain a genotype call; those below this depth are set to missing.
#' @param minQ The minimum site quality score to retain a site.
#' @param remove.completely.missing Logical. Should completely missing genotypes be removed?
#'
#' @import vcfR
#' @importFrom polyBreedR ADsplit
#'
#' @export
#'
filter_vcf <- function(x, minAD = NULL, minQ = NULL, remove.completely.missing = FALSE) {

  stopifnot(inherits(x, "vcfR"))
  merged_vcf <- x
  if (!is.null(minAD)) {
    stopifnot(is.numeric(minAD))
  } else {
    minAD <- 0
  }
  if (!is.null(minQ)) {
    stopifnot(is.numeric(minQ))
  } else {
    minQ <- 0
  }

  # Remove SNPs that are not on chromosomes
  fixed <- getFIX(merged_vcf)
  idx <- grep(pattern = "unknown", x = tolower(fixed[,"CHROM"]), ignore.case = TRUE, invert = TRUE)
  fixed <- fixed[idx, ]
  # fixed[,"CHROM"] <- str_extract(string = fixed[,"CHROM"], pattern = "chr[0-9]{1,}")
  # fixed[,"CHROM"] <- paste0("chr", str_pad(parse_number(fixed[,"CHROM"]), 2, "left", "0"))
  # fixed[,"ID"] <- paste0(fixed[,"CHROM"], "_", fixed[,"POS"])

  # filter
  merged_vcf1 <- merged_vcf[idx, ]
  merged_vcf1@fix <- cbind(fixed, INFO = getINFO(x = merged_vcf)[idx])

  # If minAD is null, do not filter
  if (!is.null(minAD)) {

    # Extract the AD information
    ad <- extract.gt(x = merged_vcf1, element = "AD")
    # Extract dp information
    dp <- extract.gt(x = merged_vcf1, element = "DP")
    # Extract the GT information; make a copy
    gt1 <- gt <- extract.gt(x = merged_vcf1, element = "GT")


    # Matrix of logicals for hets
    is_het <- gt == "0/1" | gt == "1/0"

    # Split to ref and alt
    ad_ref <- ADsplit(AD = ad, ALT = FALSE)
    ad_alt <- ADsplit(AD = ad, ALT = TRUE)

    # Apply a minimum depth filter to the AD data strings
    # First for homozy
    which_min_ad_hom <- !is_het & (ad_ref >= minAD | ad_alt >= minAD)
    which_min_ad_het <- is_het & (ad_ref >= minAD & ad_alt >= minAD)
    which_min_ad <- which_min_ad_hom | which_min_ad_het


    # If a sample is NA or does not meet the AD threshold, set as missing (./.)
    if (any(is.na(gt1))) {
      gt1[is.na(gt1)] <- "./."
    }

    if (any(!which_min_ad)) {
      gt1[!which_min_ad] <- "./."
    }

    gt_new <- matrix(paste0(gt1, ":", ad, ":", dp), nrow = nrow(gt1), ncol = ncol(gt1), dimnames = dimnames(gt1))
    gt_new <- cbind(FORMAT = "GT:AD:DP", gt_new)
    merged_vcf1@gt <- gt_new

  }

  # Filter on minQ
  if (all(is.na(merged_vcf1@fix[,"QUAL"]))) {
    idx <- seq_along(merged_vcf1@fix[,"QUAL"])
  } else {
    idx <- which(as.numeric(merged_vcf1@fix[,"QUAL"]) >= minQ)
  }
  merged_vcf1 <- merged_vcf1[idx, ]

  if (remove.completely.missing) {
    gt <- extract.gt(x = merged_vcf1, element = "GT")
    idx <- which(rowMeans(gt == "./.", na.rm = TRUE) < 1)
    merged_vcf1 <- merged_vcf1[idx, ]

  }

  return(merged_vcf1)

}



