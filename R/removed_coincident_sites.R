#' Remove SNPs with duplicate positions in a VCF
#'
#' @param x An \code{vcfR} object.
#'
#' @importFrom vcfR extract.gt
#' @importFrom polyBreedR GT2DS
#'
#' @export
#'
remove_coincident_sites <- function(x) {

  stopifnot(inherits(x, "vcfR"))

  # Pull out fix
  fix <- x@fix
  gt <- extract.gt(x = x, element = "GT", IDtoRowNames = FALSE)

  # Find SNPs with the same positions
  all_pos <- paste0(fix[,"CHROM"], fix[,"POS"])
  dup_pos <- unique(all_pos[duplicated(all_pos)])

  # If empty, return the vcf
  if (length(dup_pos) == 0) {
    cat("\nNo duplicate SNP positions found. Returning original vcfR object...")
    return(x)

  } else {
    # Index of sites to keep
    idx_keep <- which(!all_pos %in% dup_pos)
    idx_dup_keep <- numeric()

    # Iterate over duplicate sites
    for (pos in dup_pos) {
      ii <- which(all_pos %in% pos)
      gt_ii <- gt[ii,,drop = FALSE]

      # If all are the same, return the first site
      n_calls <- apply(X = gt_ii, MARGIN = 2, FUN = function(x) length(unique(x)))
      if (all(n_calls == 1)) {
        idx_dup_keep <- c(idx_dup_keep, ii[1])
      } else {
        # Find the one with the least missing data
        p_miss <- rowMeans(is.na(gt_ii))
        # If the missing rate is the same, find the one with the higher MAF
        if (length(unique(p_miss)) == 1) {
          ds <- GT2DS(GT = gt_ii, diploidize = TRUE, n.core = 1)
          maf <- calc_maf(t(ds-1), check.matrix = FALSE)
          idx_dup_keep <- c(idx_dup_keep, ii[which.max(maf)])

        } else {
          idx_dup_keep <- c(idx_dup_keep, ii[which.min(p_miss)])

        }
      }
    }

    # Total index of sites to keep
    idx_keep <- sort(c(idx_keep, idx_dup_keep))

    # Report
    cat(sub("X", replacement = nrow(fix) - length(idx_keep), x = "\nNumber of sites with duplicate genome position removed: X."))

    # Subset the VCF and return
    return(x[idx_keep, ])

  }

}
