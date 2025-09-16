#' Summarize a \code{geno} object (and derived objects)
#'
#' @param geno A \code{geno} object (or its derivatives).
#'
#' @export
#'
summarize_geno <- function(geno) {

  stopifnot(inherits(geno, "geno"))

  # Data type
  phased <- geno@phased
  phased <- ifelse(phased, "phased", "unphased")
  ploidy <- geno@ploidy
  geno_mat <- geno@geno.mat

  marNames <- colnames(geno_mat)
  indNames <- row.names(geno_mat)
  nMar <- length(marNames)
  nInd <- length(indNames)

  # Compute minor allele counts
  # First figure out which is the minor allele
  p <- colMeans(x = geno_mat, na.rm = TRUE) / ploidy # This is the frequency of the coded allele
  # If the frequency is > 0.5, then the coded allele is the major allele
  coded_is_major <- p > 0.5

  # Get counts of alleles at each locus
  geno_levels <- seq(0, ploidy)
  allele_count <- apply(geno_mat, 2, function(z) {
    z_count <- table(factor(round(z), levels = geno_levels))
    coded_count <- sum(geno_levels * z_count)
    non_coded_count <- sum((ploidy - geno_levels) * z_count)
    setNames(c(coded_count, non_coded_count), c(1, 2))
  })

  # Compute minor allele counts
  mac <- pmin(allele_count[1,], allele_count[2,])
  mac[is.na(p)] <- as.numeric(NA)

  # Compute missingness
  indMissing <- rowMeans(is.na(geno_mat))
  marMissing <- colMeans(is.na(geno_mat))

  output <- new(Class = "geno.summary",
                ploidy = as.integer(ploidy),
                nMar = nMar,
                nInd = nInd,
                mac = mac,
                mar.missing = marMissing,
                ind.missing = indMissing)

  return(output)

}
