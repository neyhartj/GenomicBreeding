#' Filter SNP marker genotypes
#'
#' @param geno A data.frame of genotype calls, where the first three columns are "marker," "chrom,", and "position."
#' Subsequent columns are genotype calls for individuals.
#' @param max.ind.miss The maximum rate of missingness to retain an individual
#' @param max.mar.miss The maximum rate of missingness to retain an marker
#' @param min.mac The minimum minor allele count to retain a marker.
#' @param max.r2 The maximum squared correlation between any pair of markers. The \code{r^2} is used directly in filtering. For each pair of markers exceeding this threshold, the one with the greater minor allele frequency is retained. Set \code{max.r2 = 1} for no filtering.
#'
#' @return
#' Variable of class \code{\link{geno}}.
#'
#' @examples
#' # Use data from the PopVar package
#' data("cranberry_geno", package = "PopVar")
#'
#' # Remove cM, rename
#' geno <- cranberry_geno[-3]
#' names(geno)[1:3] <- c("marker", "chrom", "position")
#'
#' # Filter
#' geno_in <- filter_geno
#'
#' @export
#'
filter_geno <- function(geno, max.ind.miss = 0.5, max.mar.miss = 0.5, min.mac = 5, max.r2 = 1) {

  stopifnot(is.data.frame(geno))
  cols.match <- names(geno)[1:3] == c("marker", "chrom", "position")
  if (any(!cols.match)) {
    stop("The first 3 columns of 'geno' should be 'marker', 'chrom', and 'position'.")
  }

  # Error checking
  stopifnot(max.mar.miss >= 0 & max.mar.miss <= 1)
  stopifnot(max.ind.miss >= 0 & max.ind.miss <= 1)
  stopifnot(min.mac >= 0)
  stopifnot(max.r2 >= 0 & max.r2 <= 1)

  geno_orig <- geno

  geno <- as.matrix(geno_orig[,-1:-3])
  row.names(geno) <- all.marks <- geno_orig$marker
  all.marks <- row.names(geno)
  n <- ncol(geno)
  p <- nrow(geno)
  ix <- which(apply(geno, 1, sd, na.rm = T) > 0)
  nd <- nrow(geno) - length(ix)
  cat(sub("X", nd, "Removed X markers without genetic variance\n"))
  geno <- geno[ix, , drop = FALSE]

  mar_miss <- rowMeans(is.na(geno))
  iu <- which(mar_miss <= max.mar.miss)
  nd2 <- nrow(geno) - length(iu)
  cat(sub("X", nd2, "Removed X markers due to missing data\n"))
  geno <- geno[iu, , drop = FALSE]

  ind_miss <- colMeans(is.na(geno))
  iv <- which(ind_miss <= max.ind.miss)
  ni <- ncol(geno) - length(iv)
  cat(sub("X", ni, "Removed X individuals due to missing data\n"))
  geno <- geno[, iv, drop = FALSE]

  # Calculate maf
  maf <- calc_maf(x = t(geno) - 1, check.matrix = FALSE)
  # How many markers
  iw <- which(maf >= min.mac / ncol(geno))
  nm <- nrow(geno) - length(iw)
  cat(sub("X", nm, "Removed X markers due to low minor allele count\n"))
  geno <- geno[iw, , drop = FALSE]

  # Compute pairwise r2 - only if max.r2 < 1
  if (max.r2 < 1) {
    # Impute genotypes with the mean, compute the correlation
    geno1 <- apply(X = geno, MARGIN = 1, FUN = function(snp) {
      snp[is.na(snp)] <- mean(snp, na.rm = TRUE)
      snp
    })
    geno1 <- geno1 - 1

    mar_r2 <- cor(x = geno1)^2

    # Run SNP pruning
    pruned_ans <- prune_LD(x = geno1, cor.mat = mar_r2, r2.max = max.r2, check.matrix = FALSE)
    # Report
    mars_kept <- colnames(pruned_ans)
    nr <- nrow(geno) - length(mars_kept)
    cat(sub("X", nr, "Removed X markers due to high pairwise LD\n"))

    geno <- geno[mars_kept, , drop = FALSE]

  }

  # Get the marker information back
  snp_info <- geno_orig[c(1:3)]
  snp_info <- snp_info[snp_info$marker %in% row.names(geno), ]
  geno_out <- cbind(snp_info[match(x = row.names(geno), table = snp_info$marker), ], geno)
  row.names(geno_out) <- NULL
  return(geno_out)

}
