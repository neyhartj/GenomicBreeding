#' Convert a VCF object to a geno object
#'
#' @param x A \code{vcfR} object.
#' @param ploidy The ploidy level.
#' @param n.core The number of computer cores to use.
#'
#' @importFrom vcfR extract.gt
#' @importFrom polyBreedR GT2DS
#'
#' @export
#'
vcf2genomat <- function(vcf.in, ploidy = 2, n.core = 1) {
  stopifnot(inherits(vcf.in, "vcfR"))
  stopifnot(ploidy %in% c(2, 4, 6))
  stopifnot(n.core >= 1)

  # Get GT and dosage
  gt <- extract.gt(x = vcf.in, element = "GT")
  # Detect and correct phasing
  gt <- sub("\\|", replacement = "/", x = gt)

  dos <- GT2DS(GT = gt, n.core = n.core)
  # add map information
  map_info <- as.data.frame(vcf.in@fix[,c("ID", "CHROM", "POS")])
  names(map_info) <- c("marker", "chrom", "position")
  any_na_marker <- which(is.na(map_info$marker))
  if (length(any_na_marker) > 0) {
    map_info$marker[any_na_marker] <- paste0(map_info$chrom[any_na_marker], "_", map_info$position[any_na_marker])
  }
  map_info$position <- as.numeric(map_info$position)
  map_info <- map_info[order(map_info$chrom, map_info$position), ]
  as.data.frame(cbind(map_info, dos[map_info$marker,]))
}

