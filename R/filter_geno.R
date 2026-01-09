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


#' Prune SNPs based on linkage disequilibrium
#'
#' @param x An n x p marker matrix of n individuals and p markers coded as -1, 0, and 1
#' @param cor.mat A similarity matrix between SNPs. No calculation of LD will be made if this matrix is passed.
#' for homozygous alternate, heterozygous, and homozygous reference.
#' @param r2max The maximum acceptable value of r2 (i.e. LD) or similarity between any two markers.
#' @param check.matrix Logical. Should the marker matrix 'x' be checked? Use check.matrix = FALSE
#' if 'x' contains imputed decimal genotypes.
#'
#' @return
#' A marker matrix without markers in high LD
#'
#' @export
#'
prune_LD <- function(x, cor.mat, r2.max = 0.80, check.matrix = TRUE) {

  ## Error checking
  # Check the marker matrix
  stopifnot(is.logical(check.matrix))
  if (check.matrix) stopifnot(check_marker_matrix(x))
  if (!missing(cor.mat)) {
    stopifnot(is.matrix(cor.mat))
  }
  # check r2max
  stopifnot(r2.max >= 0 & r2.max <= 1)

  # calculate minor allele frequency
  maf <- calc_maf(x, check.matrix = check.matrix)

  if (missing(cor.mat)) {

    # calculate the correlation between all markers; square it
    if (any(is.na(x))) {
      all_marker_r <- cor(x, use = "pairwise.complete.obs")^2

    } else {
      all_marker_r <- cor(x)^2
    }

  } else {
    all_marker_r <- cor.mat

  }



  # Set the lower half (including diagonal) to NA
  all_marker_r1 <- all_marker_r
  all_marker_r1[lower.tri(all_marker_r1, diag = TRUE)] <- NA

  # Get a matrix of those entries that are elevated
  elevated_marker_r <- all_marker_r1 > r2.max
  # Get the coordinates of those entries
  which_elevated <- which(x = elevated_marker_r, arr.ind = TRUE)

  markers_remove <- character()
  # While loop
  i = 1
  while(nrow(which_elevated) > 0) {

    # Subset the first row
    coords <- which_elevated[1,]
    # Extract the coordinate
    r2_coord <- all_marker_r1[coords[1], coords[2], drop = FALSE]

    # marker pair
    markers <- unlist(dimnames(r2_coord))
    # Identify the marker with higher MAF
    higher_maf <- which.max(maf[markers])

    # Toss that marker
    marker_remove <- names(higher_maf)
    markers_remove[i] <- marker_remove

    # Find the row/col containing this marker
    row_remove <- col_remove <- which(row.names(all_marker_r1) == marker_remove)

    which_elevated <- subset.matrix(which_elevated, which_elevated[,"row"] != row_remove & which_elevated[,"col"] != col_remove)

    # advance i
    i <- i + 1

  }

  # Remove the markers from the marker matrix
  cols_keep <- setdiff(seq_len(ncol(x)), which(colnames(x) %in% markers_remove))
  x[,cols_keep,drop = FALSE]

}
