#' Summarize a \code{geno} object (and derived objects)
#'
#' @export
#'
setMethod(
  "print",
  signature(x = "geno.summary"),
  function(x) {
    cat("Summary of a geno object:\n")

    cat("\nNumber of individuals:", x@nInd)
    cat("\nNumber of markers:", x@nMar, "\n")

    # Compute minor allele frequency
    maf <- x@mac / (x@nInd * x@ploidy)

    # Quantiles of MAF
    maf_quants <- round(quantile(maf, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE), 3)
    names(maf_quants) <- c("Min", "1st Qu.", "Median", "3rd Qu.", "Max")
    cat("\nMinor allele frequency (MAF) summary:\n")
    print(maf_quants)

    # Quantiles of marker missingness
    mar_miss_quants <- round(quantile(x@mar.missing, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE), 3)
    names(mar_miss_quants) <- c("Min", "1st Qu.", "Median", "3rd Qu.", "Max")
    cat("\nMarker missingness proportion summary:\n")
    print(mar_miss_quants)

    # Quantiles of individual missingness
    ind_miss_quants <- round(quantile(x@ind.missing, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE), 3)
    names(ind_miss_quants) <- c("Min", "1st Qu.", "Median", "3rd Qu.", "Max")
    cat("\nIndividual missingness proportion summary:\n")
    print(ind_miss_quants)

  }
)


#' @describeIn summarize_geno
#'
#' @importFrom graphics plot
#'
#' @export
#'
setMethod(
  "plot",
  signature(x = "geno.summary"),
  function(x) {

    # Save current plot params
    old_par <- par(no.readonly = TRUE)
    # Set plot params for 2x2 layout
    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

    # Plot MAF distribution
    maf <- x@mac / (x@nInd * x@ploidy)
    hist(maf, main = "Minor Allele Frequency (MAF) Distribution", xlab = "MAF", breaks = 20, col = "lightblue", border = "black")

    # Plot marker missingness distribution
    hist(x@mar.missing, main = "Marker Missingness Distribution", xlab = "Proportion Missing", breaks = 20, col = "lightgreen", border = "black")

    # Plot individual missingness distribution
    hist(x@ind.missing, main = "Individual Missingness Distribution", xlab = "Proportion Missing", breaks = 20, col = "lightcoral", border = "black")

    # Revert the plot params
    par(old_par)

  }
)



#' Print a geno object
#'
#' @export
#'
setMethod(
  "print",
  signature(x = "geno"),
  function(x) {
    class <- class(x)
    cat("This is an object of class", class, "\n")

    phased <- x@phased
    phased <- ifelse(phased, "phased", "unphased")
    ploidy <- x@ploidy
    geno_mat <- x@geno.mat
    nInd <- nrow(geno_mat)
    nMar <- ncol(geno_mat)

    cat("\nMarker data is", phased)
    cat("\nPloidy level is", ploidy, "\n")

    cat("\nNumber of individuals =", nInd)
    cat("\nNumber of markers =", nMar, "\n")

    cat("\nUse 'summarize_geno()' to get more detailed information about this object.\n")

  }
)
