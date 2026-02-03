#' Convert rrblup formatted genotype data to a geno matrix
#'
#' @param x A data.frame of genotype calls, where the first three columns are "marker," "chrom,", and "position."
#' Subsequent columns are genotype calls for individuals.
#' @param transpose Logical. Should the output matrix be indiv x markers (TRUE) or markers x indiv (FALSE)
#' @param output The output type, either 0-1-2 for homozygous reference, heterozygous, homozygous alternate, or -1-0-1 for the same.
#'
#' @export
#'
rrblup2genomat <- function(x, transpose = TRUE, output = c("012", "-101")) {
  stopifnot(is.logical(transpose))
  output <- match.arg(output)
  x <- as.data.frame(x)
  x1 <- x[,-1:-3]
  row.names(x1) <- x[[1]]
  x2 <- as.matrix(x1)
  elements <- na.omit(unique(as.vector(x2)))
  if (all(elements %in% c(0, 1, 2))) {
    input <- "012"
  } else if (all(elements %in% c(-1, 0, 1))) {
    input <- "-101"
  } else {
    stop("Input format not recognized")
  }
  if (input == output) {
    x2 <- x2
  } else if (input == "012" & output == "-101") {
    x2 <- x2 - 1
  } else if (input == "-101" & output == "012") {
    x2 <- x2 + 1
  }

  if (transpose) {
    return(t(x2))
  } else {
    return(x2)
  }
}
