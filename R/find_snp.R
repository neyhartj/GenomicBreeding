#' Find the SNP position between two sequences
#'
#' @description
#' This is a naive function to find nucleotide mismatches between two sequences.
#'
#' @param seq1 A sequence
#' @param seq2 Another sequence
#' @param unambiguous Logical. Use only unambiguous nucleotides.
#'
#'
#' @export
#'
find_snp <- function(seq1, seq2, unambiguous = FALSE) {
  # Ensure both sequences are the same length
  if (nchar(seq1) != nchar(seq2)) {
    stop("Sequences must be of equal length.")
  }

  # Split sequences into individual characters
  seq1_split <- strsplit(seq1, "")[[1]]
  seq2_split <- strsplit(seq2, "")[[1]]

  # Identify positions where the sequences differ
  snp_positions <- which(seq1_split != seq2_split)

  if (unambiguous) {
    snp_positions_ref <-  snp_positions[which(sapply(X = snp_positions, FUN = function(x) seq1_split[x]) %in% c("A", "C", "G", "T"))]
    snp_positions_alt <-  snp_positions[which(sapply(X = snp_positions, FUN = function(x) seq2_split[x]) %in% c("A", "C", "G", "T"))]
    snp_positions <- intersect(snp_positions_ref, snp_positions_alt)
  }

  # Return the SNP positions
  if (length(snp_positions) == 0) {
    return(numeric(0))
  } else {
    return(snp_positions)
  }
}
