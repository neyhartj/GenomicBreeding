#' Load the gene models for the cranberry Ben Lear reference genome
#'
#' @param path A path to the GFF file (or url). If NULL (default), the raw GFF will be downloaded from \link{vaccinium.org}.
#' @param download Should the GFF file be downloaded?
#'
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges end
#' @importFrom Seqinfo seqnames `seqlengths<-`
#'
#' @export
#'
load_cran_gff <- function(path = NULL, download = FALSE) {

  if (is.null(path)) {
    path <- "https://www.vaccinium.org/vaccinium_downloads/Vaccinium_macrocarpon/Vmacrocarpon_BenLear_genome_v1/annotation/Vmacrocarpon_BenLear_v1-annotation.gff3"
    if (!download) stop("If 'path' is NULL, set 'download' to TRUE.")
    gff_data <- import(path, genome = "BenLear")

  } else {
    gff_data <- import(path)
  }

  chr_rang <- range(gff_data)
  maxlens <- tapply(X = end(chr_rang), INDEX = seqnames(chr_rang), FUN = max)
  maxlens <- setNames(as.numeric(maxlens), names(maxlens))
  seqlengths(gff_data) <- maxlens

  # Return the db
  return(gff_data)
}



#' Find genes from a table of markers
#'
#' @description
#' Find genes based on a table of markers, such as significant QTL peaks from a mapping
#' study or GWAS, or outlier SNPs from other analyses.
#'
#' @param qtl A data frame with the name, chromosome, and physical position of markers; one row per marker.
#' @param genes A \code{GRanges} object with gene model data for the same reference genome as the position of the markers in \code{qtl}.
#' @param window Numeric size of the window from each marker in \code{qtl} to find genes.
#'
#'
#' @importFrom Seqinfo seqnames seqlengths `seqlengths<-`
#' @importFrom IRanges IRanges findOverlaps
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors queryHits subjectHits mcols
#'
#' @export
#'
find_genes <- function(qtl, genes, window = 2e4) {
  # Error checking
  stopifnot(is.data.frame(qtl))
  if (ncol(qtl) < 3) stop ("'qt' should be a data frame with 3 columns.")
  stopifnot(inherits(genes, "GRanges"))
  stopifnot(is.numeric(window))
  stopifnot(window >= 0)

  # Create a GRanges object from qtl
  qtl2 <- GRanges(
    seqnames = qtl[[2]],
    ranges = IRanges(start = qtl[[3]], end = qtl[[3]]),
    marker = qtl[[1]]
  )

  seqlengths(qtl2) <- seqlengths(genes)[names(seqlengths(qtl2))]

  # Get genes only from the genes input -- CHANGE THIS LATER
  all_genes <- genes[genes$type == "gene", ]
  all_trans <- genes[genes$type == "mRNA", ]

  # Find overlaps
  overlaps <- findOverlaps(query = qtl2, subject = all_genes, maxgap = window)

  results_df <- data.frame(
    marker   = qtl2$marker[queryHits(overlaps)],
    chrom    = seqnames(qtl2)[queryHits(overlaps)],
    position    = start(qtl2)[queryHits(overlaps)],

    gene_id    = all_genes$ID[subjectHits(overlaps)],
    gene_chr   = seqnames(all_genes)[subjectHits(overlaps)],
    gene_start = start(all_genes)[subjectHits(overlaps)],
    gene_end   = end(all_genes)[subjectHits(overlaps)]
  )

  results_df$distance_to_start <- abs(results_df$position - results_df$gene_start)
  results_df$distance_to_end <- abs(results_df$position - results_df$gene_end)

  # Merge annotation data in
  anno_df <- as.data.frame(mcols(all_trans))

  results_df1 <- merge(x = results_df, y = anno_df, by.x = "gene_id", by.y = "Parent")

  # Return this:
  return(results_df1)

}






















