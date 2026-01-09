#' Filter MADC data
#'
#' @description
#' Applies filters to MADC data from DArTag experiments
#'
#' @param madc A data.frame of MADC data (read from a .csv); ignored if \code{madc.file} is passed.
#' @param madc.file The name of the MADC file; ignored if \code{madc} is passed.
#' @param min.mh.reads Minimum number of reads at genotyped samples to retain a RefMatch or AltMatch allele.
#' @param min.mh.samples Minimum number of samples required to retain a RefMatch or AltMatch allele.
#' has at least \code{min.mh.reads} reads.
#' @param min.mh.ref.identity Minimum proportion of identity of a RefMatch or AltMatch allele to the reference.
#' @param min.mh.ref.coverage Minimum proportion of coverage of a RefMatch or AltMatch allele to the reference.
#' @param min.locus.samples Minimum number of samples to retain a locus.
#' @param min.locus.reads Minimum number of reads at a genotyped locus to retain a locus.
#' @param max.locus.haps Maximum number of haplotypes at a locus. If greater, only the Ref and Alt haplotypes are retained.
#' @param max.sample.missing Maximum missing data proportion to retain a sample/individual.
#'
#' @import polyRAD
#'
#' @export
#'
filter_madc <- function(madc, madc.file, min.mh.reads = 2, min.mh.samples = 10, min.mh.ref.identity = 0.90,
                        min.mh.ref.coverage = 0.90, min.locus.samples = 10, min.locus.reads = 5,
                        max.locus.haps = 10, max.sample.missing = 0.95) {

  if (!missing(madc)) {
    stopifnot(is.data.frame(madc))
    madc <- as.data.frame(madc)
  } else if (!missing(madc.file)) {
    stopifnot(file.exists(madc.file))
    madc <- read.csv(file = madc.file, as.is = TRUE, check.names = FALSE)
  } else {
    stop("One of 'madc' or 'madc.file' must be passed.")
  }

  stopifnot(all(c("AlleleID", "CloneID", "AlleleSequence") %in% names(madc)))

  # Summarize and print
  n_loci <- length(unique(madc$CloneID))
  n_alleles <- nrow(madc)
  n_samples <- ncol(madc) - 3
  n_alleles_per_loci <- sapply(split(madc$AlleleID, madc$CloneID), length)

  cat("\nSummary of the input MADC file:")
  cat("\nNumber of samples:", n_samples)
  cat("\nNumber of loci:", n_loci)
  cat("\nTotal number of microhaplotype alleles:", n_alleles)
  cat("\nDistribution of microhaplotype alleles per locus:\n")
  print(summary(n_alleles_per_loci))

  cat("\n\nFiltering the MADC file...")

  # Create a matrix of allele read counts
  madc_mat <- as.matrix(madc[,-1:-3])
  row.names(madc_mat) <- madc$AlleleID

  # First create an index of the ref and alt haplotypes; these will always be retained
  idx_refalt <- grep(pattern = "Ref_[0]{1,}1|Alt_[0]{1,}2", x = row.names(madc_mat), value = TRUE)

  # Compute the number of samples per haplotype that are above the read count
  hap_sample_read_ct <- rowSums(madc_mat >= min.mh.reads)
  # Determine which haplotypes meet the criteria
  which_hap_sample_read_ct <- which(hap_sample_read_ct >= min.mh.samples)

  # Split madc sequence information by locus
  madc_by_locus <- split(madc[,1:3], madc$CloneID)

  # Calculate the percent identity of the haplotype sequences to the reference
  seqs_by_locus <- lapply(X = madc_by_locus, function(x) setNames(x$AlleleSequence, x$AlleleID) )
  perc_ident_seqs <- map(seqs_by_locus, ~{
    ref_idx <- grep(pattern = "Ref_[0]{1,}1", x = names(.x))
    alt_idx <- grep(pattern = "Alt_[0]{1,}2", x = names(.x))
    match_idx <- grep(pattern = "Match", x = names(.x))
    other_idx <- setdiff(seq_along(.x), c(ref_idx, alt_idx, match_idx))
    dist <- polyRAD:::.nucdist(alleles1 = .x[c(match_idx, other_idx)], alleles2 = .x[ref_idx])
    len_ref <- nchar(.x[ref_idx])
    perc_ident <- 1 - (dist / len_ref)
    # Return
    matrix(c(1, 1, perc_ident), ncol = 1, dimnames = list(names(.x)[c(ref_idx, alt_idx, match_idx, other_idx)], "perc_ident"))
  })

  # Determine which haplotypes meet the criteria
  perc_ident_mat <- do.call("rbind", perc_ident_seqs)
  which_hap_perc_ident <- which(perc_ident_mat[,1] >= min.mh.ref.identity)

  # apply the first round of filters
  haps_keep1 <- sort(union(idx_refalt, intersect(names(which_hap_sample_read_ct), names(which_hap_perc_ident))))
  madc1 <- madc[,1:3]
  madc1 <- madc1[madc1$AlleleID %in% haps_keep1, ]

  # Split again by locus
  madc1_by_locus <- split(madc1[,1:3], madc1$CloneID)
  # Count the number of microhaplotypes by locus; if excessive; only retain the ref/alt alleles
  haps_keep2 <- lapply(madc1_by_locus, function(x) {
    ref_alt_idx <- grep(pattern = "Ref_[0]{1,}1|Alt_[0]{1,}2", x = x$AlleleID)
    n_hap_alleles <- nrow(x) - length(ref_alt_idx)
    if (n_hap_alleles > max.locus.haps) {
      x$AlleleID[ref_alt_idx]
    } else {
      x$AlleleID
    }
  })
  haps_keep2 <- unlist(haps_keep2)

  # apply the second round of filters
  madc2 <- madc1[madc1$AlleleID %in% haps_keep2, ]

  # Create the full matrix
  madc_mat1 <- madc_mat[madc2$AlleleID, , drop = FALSE]
  # Split alleles by locus
  alleles_by_locus <- split(madc2$AlleleID, madc2$CloneID)
  # For each locus, compute the number of samples with the minimum total read count
  locus_sample_read_ct <- sapply(alleles_by_locus, FUN = function(x) {
    dat <- madc_mat1[x, , drop = FALSE]
    sum(colSums(dat) >= min.locus.reads)
  })

  # Determine which loci meet the criteria
  which_locus_sample_read_ct <- which(locus_sample_read_ct >= min.locus.samples)

  # Apply the third round of filters
  loci_keep1 <- which(sapply(strsplit(row.names(madc_mat1), "\\|"), "[[", 1) %in% names(which_locus_sample_read_ct))
  madc_mat2 <- madc_mat1[loci_keep1, , drop = FALSE]
  madc3 <- madc2[madc2$AlleleID %in% row.names(madc_mat2), ]

  # Split alleles by locus
  alleles_by_locus <- split(madc3$AlleleID, madc3$CloneID)
  # For each locus, create a matrix of TRUE/FALSE whether a genotype has any data at the locus
  marker_genotype_rate_mat <- lapply(alleles_by_locus, FUN = function(x) {
    dat <- madc_mat2[x, , drop = FALSE]
    tot_reads_by_idv <- colSums(dat)
    tot_reads_by_idv > 0
  })
  marker_genotype_rate_mat <- do.call("rbind", marker_genotype_rate_mat)

  # Determine the missing rate per genotype
  sample_missing_rate <- 1 - colMeans(marker_genotype_rate_mat)
  # Remove samples with excess missing data
  samples_keep1 <- which(sample_missing_rate <= max.sample.missing)

  madc_mat3 <- madc_mat2[, names(samples_keep1), drop = FALSE]

  # Create a new madc data.frame
  madc_out <- cbind(madc3, madc_mat3[madc3$AlleleID, , drop = FALSE])

  # Summarize the filtered madc
  # Summarize and print
  n_loci <- length(unique(madc_out$CloneID))
  n_alleles <- nrow(madc_out)
  n_samples <- ncol(madc_out) - 3
  n_alleles_per_loci <- sapply(split(madc_out$AlleleID, madc_out$CloneID), length)

  cat("\n\nSummary of the filtered MADC file:")
  cat("\nNumber of samples:", n_samples)
  cat("\nNumber of loci:", n_loci)
  cat("\nTotal number of microhaplotype alleles:", n_alleles)
  cat("\nDistribution of microhaplotype alleles per locus:\n")
  print(summary(n_alleles_per_loci))

  return(madc_out)

}
