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



