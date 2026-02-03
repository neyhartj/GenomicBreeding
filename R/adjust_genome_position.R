#' Adjust the genome position of target and off-target DArT markers
#'
#' @description
#' Use BLAST alignment results to modify the genome position of DArT SNP markers.
#'
#' @param vcf.in A \code{vcfR} object with SNP markers (i.e. from the output of \code{\link{BIGr}[madc2vcf_all]}).
#' @param alignment A data frame with information for genome position conversion.
#' @param reference.name The name of the reference genome to which positions will be converted.
#' @param locus.col The name of the column in \code{alignment} that contains the target locus ID.
#' @param target.snp.pos.col The name of the column in \code{alignment} that contains the position of the target SNP in the DArT tag.
#' @param genome.chrom.col The name of the column in \code{alignment} that contains the new genome chromosome names.
#' @param locus.start.col The name of the column in \code{alignment} that contains the start position of the query tag.
#' @param locus.end.col The name of the column in \code{alignment} that contains the end position of the query tag.
#' @param genome.start.col The name of the column in \code{alignment} that contains the start position in the genome of the alignment.
#' @param genome.end.col The name of the column in \code{alignment} that contains the end position in the genome of the alignment.
#'
#' @details
#' The function will remove sites that do not have a corresponding position in the new genome.
#'
#'
#' @return
#' A \code{vcfR} object with SNP positions converted to the position within the genome of interest.
#'
#' @importFrom vcfR extract.info
#'
#' @export
#'
#'
adjust_genome_position <- function(vcf.in, alignment, reference.name, locus.col = "locus_id", target.snp.pos.col = "snp_position_in_locus",
                                   genome.chrom.col = "benlear_v1_chrom", locus.start.col = "query_benlear_v1_start", locus.end.col = "query_benlear_v1_end",
                                   genome.start.col = "benlear_v1_start", genome.end.col = "benlear_v1_end") {

  # Error handling
  stopifnot(inherits(vcf.in, "vcfR"))
  stopifnot(is.data.frame(alignment))
  alignment <- as.data.frame(alignment)
  stopifnot(is.character(reference.name))
  # Check columns in data frame
  col_names <- c("locus.col" = locus.col, "target.snp.pos.col" = target.snp.pos.col, "genome.chrom.col" = genome.chrom.col, "locus.start.col" = locus.start.col,
                 "locus.end.col" = locus.end.col, "genome.start.col" = genome.start.col, "genome.end.col" = genome.end.col)
  col_missing <- !col_names %in% names(alignment)
  if (any(col_missing)) stop("The following provided column names were not found in the 'alignment' input:", paste0(col_names[col_missing], collapse = ", "))

  # Select the relevant columns from alignment
  alignment1 <- alignment[col_names]

  ## Change positions
  # Get the SNP IDs from the vcf
  fix_in <- vcf.in@fix
  marker_ids <- fix_in[, "ID"]
  # Make sure target information is present in fix_in
  target_data <- extract.info(x = vcf.in, element = "TYPE")
  if (all(is.na(target_data))) stop("Target vs off-target information (coded as TYPE in INFO) is not present in 'vcf.in'.")

  # Get the target SNPs
  target_snps_idx <- which(target_data == "target")
  target_snps <- marker_ids[target_snps_idx]

  # Empty framework for the new fix object
  new_fix <- data.frame(orig_id = marker_ids, orig_chrom = as.character(NA), orig_pos = as.numeric(NA),
                        target = target_data == "target", dist.to.closest.snp = as.numeric(NA), closest.target.snp = as.character(NA))
  marker_ids_split <- strsplit(x = marker_ids, split = "_")
  new_fix$orig_chrom <- sapply(marker_ids_split, "[[", 1)
  new_fix$orig_pos <- as.numeric(sapply(marker_ids_split, "[[", 2))

  # Create a similar file for the alignment data
  alignment_fix <- alignment1[locus.col]
  alignment_fix_split <- strsplit(x = alignment_fix[[1]], split = "_")
  alignment_fix$assumed.chrom <- sapply(alignment_fix_split, "[[", 1)
  alignment_fix$assumed.pos <- as.numeric(sapply(alignment_fix_split, "[[", 2))

  # Iterate over rows in new_fix
  for (i in seq_len(nrow(new_fix))) {
    # subset rows in alignment_fix from the same chromosome
    alignment_fix_chrom_i <- alignment_fix[alignment_fix$assumed.chrom == new_fix$orig_chrom[i], ]
    dist_i <- new_fix$orig_pos[i] - alignment_fix_chrom_i$assumed.pos
    dist_i_choose <- which.min(abs(dist_i))

    # Find the closest target SNP
    new_fix$dist.to.closest.snp[i] <- dist_i[dist_i_choose]
    new_fix$closest.target.snp[i] <- alignment_fix_chrom_i$locus_id[dist_i_choose]
  }

  # Merge with alignment fix
  new_fix1 <- merge(x = new_fix, y = alignment1, by.x = "closest.target.snp", by.y = locus.col)
  # Determine the reverse strand alignments
  new_fix1$is.reverse <- new_fix1[[genome.start.col]] > new_fix1[[genome.end.col]]
  # Determine SNP position in the genome
  genome_snp_pos_forward <- new_fix1[[genome.start.col]] + new_fix1[[target.snp.pos.col]] - new_fix1[[locus.start.col]] + new_fix1$dist.to.closest.snp
  genome_snp_pos_reverse <- new_fix1[[genome.start.col]] - new_fix1[[target.snp.pos.col]] + new_fix1[[locus.start.col]] - new_fix1$dist.to.closest.snp

  new_fix1$genome.snp.pos <- ifelse(new_fix1$is.reverse, genome_snp_pos_reverse, genome_snp_pos_forward)

  # Index of loci to keep
  idx_keep <- which(!is.na(new_fix1[[genome.chrom.col]]))

  # Create a new fix with the new SNP positions
  fix1 <- cbind(CHROM = new_fix1[[genome.chrom.col]], POS = as.character(new_fix1$genome.snp.pos),
                ID = paste0(new_fix1[[genome.chrom.col]], "_", str_pad(string = as.character(new_fix1$genome.snp.pos), width = 9, side = "left", pad = "0")),
                REF = fix_in[,"REF"], ALT = fix_in[,"ALT"], QUAL = fix_in[,"QUAL"], FILTER = fix_in[,"FILTER"],
                # Add the original SNP ID to the INFO string
                INFO = paste0(fix_in[,"INFO"], ";ORIGID=", fix_in[,"ID"]))

  fix1 <- fix1[idx_keep, ]

  # Subset GT
  gt1 <- vcf.in@gt[idx_keep, ]

  # Adjust meta
  meta1 <- c("##fileformat=VCFv4.3",
             paste0("##filedate=", format(Sys.Date(), "%Y%m%d")),
             paste0("##reference=", reference.name),
             paste0("##contig=<ID=", unique(fix1[,"CHROM"]), ">"),
             grep(pattern = "INFO", x = vcf.in@meta, value = TRUE),
             "##INFO=<ID=TYPE,Type=String,Description=\"Whether the SNP is a target or off-target marker\">",
             "##INFO=<ID=ORIGID,Type=String,Description=\"The original ID of the marker\">",
             vcf.in@meta[-seq(max(grep(pattern = "INFO", x = vcf.in@meta)))]
  )

  # Make the new VCF
  vcf.out <- vcf.in
  vcf.out@meta <- meta1
  vcf.out@fix <- fix1
  vcf.out@gt <- gt1

  return(vcf.out)

}
