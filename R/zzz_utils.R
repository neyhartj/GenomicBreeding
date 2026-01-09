#' Hidden functions
#'
check_marker_matrix <- function(x) {
  # First check if it's a matrix
  if (!is.matrix(x)) {
    return(FALSE)

  } else if (!all(na.omit(unique(as.vector(x))) %in% c(-1, 0, 1), na.rm = TRUE)) {
    return(FALSE)

  } else if (any(sapply(dimnames(x), is.null))) {
    return(FALSE)

  } else {
    return(TRUE)

  }

}


#' Calculate per-snp or per-genotype missingness
#'
#' @param x An n x p marker matrix of n individuals and p markers coded as -1, 0, and 1
#' for homozygous alternate, heterozygous, and homozygous reference.
#' @param type Calculate individual or marker missingness
#' @param check.matrix Logical. Should the marker matrix 'x' be checked? Use check.matrix = FALSE
#' if 'x' contains imputed decimal genotypes.
#'
#' @export
#'
calc_missing <- function(x, type = c("individual", "marker"), check.matrix = TRUE) {

  ## Error checking
  stopifnot(is.logical(check.matrix))
  if (check.matrix) stopifnot(check_marker_matrix(x))

  miss <- is.na(x)

  # Calculate per-snp missingness
  snp_miss <- colMeans(miss)
  # Genotype missingness
  geno_miss <- rowMeans(miss)

  if (type == "individual") {
    return(geno_miss)
  } else {
    return(snp_miss)
  }

}


#' Calculate minor allele frequency
#'
#' @param x An n x p marker matrix of n individuals and p markers coded as -1, 0, and 1
#' for homozygous alternate, heterozygous, and homozygous reference.
#' @param check.matrix Logical. Should the marker matrix 'x' be checked? Use check.matrix = FALSE
#' if 'x' contains imputed decimal genotypes.
#'
#' @export
#'
calc_maf <- function(x, check.matrix = TRUE) {

  ## Error checking
  stopifnot(is.logical(check.matrix))
  if (check.matrix) stopifnot(check_marker_matrix(x))

  # Filter for MAF
  af <- colMeans(x + 1, na.rm = TRUE) / 2
  maf <- pmin(af, 1 - af)

  return(maf)

}
