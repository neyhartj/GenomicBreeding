#' Convert a rrBLUP formatted object to AlphaGenes format
#'
#' @param geno.df A data.frame of genotype calls, where the first three columns are "marker," "chrom,", and "position."
#' Subsequent columns are genotype calls for individuals.
#' @param outdir The path the output directory to save the files.
#' @param sep.by.chrom Logical. Should markers be separate by chromosome?
#' @param prefix The filename prefix
#'
#' @importFrom vcfR extract.gt getFIX
#' @importFrom polyBreedR GT2DS
#'
#' @export
#'
rrblup2ag <- function(geno.df, outdir = ".", sep.by.chrom = TRUE, prefix = "cran_snps_alphaplantimpute_input") {

  stopifnot(inherits(geno.df, "data.frame"))
  stopifnot(dir.exists(outdir))
  stopifnot(is.logical(sep.by.chrom))

  info <- geno.df[,c("marker", "chrom", "position")]
  gt <- geno.df[,setdiff(colnames(geno.df), c("marker", "chrom", "position"))]
  # Get a list of marker IDs per chromosome
  ids_per_chrom <- split(info$marker, info$chrom)
  geno1 <- t(gt)
  colnames(geno1) <- info$marker
  geno1[is.na(geno1)] <- as.numeric(9)

  outfiles <- c()

  # Split by chromosome, if called for
  if (sep.by.chrom) {
    for (i in seq_along(ids_per_chrom)) {
      chr_name <- names(ids_per_chrom)[i]
      markers <- ids_per_chrom[[i]]
      genos_i <- geno1[,markers, drop = FALSE]
      genos_i <- as.data.frame(genos_i)
      genos_i <- cbind(geno_id = row.names(genos_i), genos_i)
      outfile <- file.path(outdir, paste0(prefix, "_", chr_name, ".geno"))
      outfiles[i] <- outfile
      write.table(x = genos_i, file = outfile, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

    }

  } else {
    geno2 <- as.data.frame(geno1)
    geno2 <- cbind(geno_id = row.names(geno2), geno2)
    outfile <- file.path(outdir, paste0(prefix, ".geno"))
    outfiles <- outfile
    write.table(x = geno2, file = outfile, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

  }

  return(outfiles)

}
