#' Filter SNP marker genotypes
#'
#' @param geno A object of \class{geno} (i.e. output of function \code{\link{read_geno}}).
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
#' data("barley_geno", package = "PopVar")
#' data("cranberry_geno", package = "PopVar")
#' data("cranberry_geno_phased", package = "PopVar")
#'
#' # Read in marker genotype data
#' geno_in <- read_geno(geno = barley_geno)
#'
#' # Filter
#' geno_in <- filter_geno
#'
#' @export
#'
filter_geno <- function(geno, max.ind.miss = 0.5, max.mar.miss = 0.5, min.mac = 5, max.r2 = 1) {

  # Error checking
  stopifnot(inherits(geno_in, "geno"))
  stopifnot(max.mar.miss >= 0 & max.mar.miss <= 1)
  stopifnot(max.ind.miss >= 0 & max.ind.miss <= 1)
  stopifnot(min.mac >= 0)
  stopifnot(max.r2 >= 0 & max.r2 <= 1)

  geno_mat <- geno@geno.mat
  map <- geno@map
  haplo_mat <- geno@haplo.mat
  ploidy <- geno@ploidy

  # Run summarize_geno to get missingness and mac
  geno_sum <- summarize_geno(geno)

  ## Create base indices of markers and individuals to retain
  # index of individuals
  idx <- seq_len(nrow(geno_mat))
  # Index of markers
  idy <- seq_len(ncol(geno_mat))

  # Filter individuals by missingness
  ind_missing <- geno_sum@ind.missing
  idx_remove <- which(ind_missing > max.ind.miss)
  idx <- setdiff(idx, idx_remove)
  # Print a message
  if (length(idx) == 0) {
    stop("No individuals passed the 'max.ind.miss' threshold.")
  } else {
    cat(sub("X", length(idx), "Number of individuals passing the 'max.ind.miss' threshold = X\n"))
  }

  mar_missing <- geno_sum@mar.missing
  idy_remove <- which(mar_missing > max.mar.miss)
  idy <- setdiff(idy, idy_remove)
  # Print a message
  if (length(idy) == 0) {
    stop("No markers passed the 'max.mar.miss' threshold.")
  } else {
    cat(sub("X", length(idy), "Number of markers passing the 'max.mar.miss' threshold = X\n"))
  }

  mac <- geno_sum@mac
  idy_remove <- which(mac < min.mac)
  idy <- setdiff(idy, idy_remove)
  # Print a message
  if (length(idy) == 0) {
    stop("No markers passed the 'min.mac' threshold.")
  } else {
    cat(sub("X", length(idy), "Number of markers passing the 'min.mac' threshold = X\n"))
  }

  # Compute pairwise r2 - only if max.r2 < 1
  if (max.r2 < 1) {
    # subset the geno mat
    geno_mat1 <- geno_mat[idx, idy, drop = FALSE]

    mar_r2 <- cor(x = geno_mat1, use = "pairwise.complete.obs")^2
    # Create a matrix from this
    mar_r2[lower.tri(mar_r2, diag = TRUE)] <- NA

    mar_r2_long <- which(!is.na(mar_r2), arr.ind = TRUE)
    mar_r2_long <- cbind(mar_r2_long, r2 = mar_r2[mar_r2_long])

    # While loop to remove markers until all pairwise r2 <= max.r2
    idy_remove <- numeric(0)
    while(any(mar_r2_long[,3] > max.r2)) {
      # Find the pair with the highest correlation
      pair_test <- mar_r2_long[which.max(mar_r2_long[,3]), , drop = FALSE]
      # Find the marker with the highest MAC
      mar_test <- c(row.names(mar_r2)[pair_test[,1]], colnames(mar_r2)[pair_test[,2]])
      remove <- names(which.min(mac[mar_test]))

      # Find the position of this marker in the larger matrix
      remove_j <- which(colnames(geno_mat) == remove)
      idy_remove <- c(idy_remove, remove_j)

      # Remove that marker from the long correlation matrix
      mar_r2_long <- mar_r2_long[!(mar_r2_long[,1] == pair_test[,1] | mar_r2_long[,2] == pair_test[,2]), , drop = FALSE]

    }

    # Add the remove markers to the list
    idy <- setdiff(idy, idy_remove)

    # Print a message
    if (length(idy) == 0) {
      stop("No markers passed the 'max.r2' threshold.")
    } else {
      cat(sub("X", length(idy), "Number of markers passing the 'max.r2' threshold = X\n"))
    }

  }



  output <- new(Class = "PopVar.geno", ploidy = as.integer(ploidy), map = map, geno.mat = geno_mat, haplo.mat = haplo_mat,
                phased = phased, G = G, D = D)

  return(output)
}
