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
