#' Convert a VCF object to AlphaGenes format
#'
#' @param vcf A \code{vcfR} object.
#' @param outdir The path the output directory to save the files.
#' @param sep.by.chrom Logical. Should markers be separate by chromosome?
#' @param n.core The number of computer cores to use.
#'
#' @importFrom vcfR extract.gt getFIX
#' @importFrom polyBreedR GT2DS
#'
#' @export
#'
vcf2ag <- function(vcf, outdir = ".", sep.by.chrom = TRUE, n.core = 1) {

  stopifnot(inherits(vcf, "vcfR"))
  stopifnot(dir.exists(outdir))
  stopifnot(is.logical(sep.by.chrom))
  # Extract the GT matrix
  gt <- extract.gt(x = vcf, element = "GT")

  # Extract the info
  info <- as.data.frame(getFIX(vcf))
  # Get a list of marker IDs per chromosome
  ids_per_chrom <- split(info$ID, info$CHROM)

  # Detect if there is phasing
  if (any(grepl(pattern = "/", gt))) {
    gt <- sub("\\|", replacement = "/", x = gt)
    gt <- GT2DS(GT = gt, n.core = n.core)
    geno1 <- t(gt)
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
        outfile <- file.path(outdir, paste0("alphaplantimpute_input_", chr_name, ".geno"))
        outfiles[i] <- outfile
        write.table(x = genos_i, file = outfile, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

      }

    } else {
      geno2 <- as.data.frame(geno1)
      geno2 <- cbind(geno_id = row.names(geno2), geno2)
      outfile <- file.path(outdir, "alphaplantimpute_input.geno")
      outfiles <- outfile
      write.table(x = geno2, file = outfile, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

    }

  } else {
    # Else create library files
    marks <- vcf@fix[,"ID"]
    inds <- colnames(gt)
    # Split haplotypes
    haps <- apply(X = gt, MARGIN = 2, FUN = function(x) {
      x1 <- strsplit(x = x, split = "\\|")
      hap1 <- sapply(X = x1, FUN = "[[", 1)
      hap2 <- sapply(X = x1, FUN = "[[", 2)
      list(rbind(hap1, hap2))
    })
    haps <- lapply(X = haps, "[[", 1)
    ref_df <- mapply(haps, inds, FUN = function(x, y) `rownames<-`(x, rep(y, nrow(x))), SIMPLIFY = FALSE)
    ref_df <- do.call(rbind, ref_df)
    ref_df <- as.data.frame(cbind(id = row.names(ref_df), ref_df))

    outfiles <- c()

    # Split by chromosome, if called for
    if (sep.by.chrom) {
      for (i in seq_along(ids_per_chrom)) {
        chr_name <- names(ids_per_chrom)[i]
        markers <- ids_per_chrom[[i]]
        genos_i <- ref_df[c("id", markers)]
        outfile <- file.path(outdir, paste0("alphaplantimpute_haplotype_library_", chr_name, ".phase"))
        outfiles[i] <- outfile
        write.table(x = genos_i, file = outfile, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

      }

    } else {
      genos_i <- ref_df[c("id", markers)]
      outfile <- file.path(outdir, "alphaplantimpute_haplotype_library.phase")
      outfiles <- outfile
      write.table(x = genos_i, file = outfile, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

    }

  }

  return(outfiles)

}
