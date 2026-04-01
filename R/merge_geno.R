#' Merge SNP marker genotypes
#'
#' @param geno.list A list of data frames of genotype calls, where the first three columns are "marker," "chrom,", and "position."
#' @param run.ids A vector of characters defining the run ID for each data frame in \code{geno.list}.
#' This ID will be added as a suffix in repeated sample ID in case they exist in different files.
#' @param sep The string to separate sample names from the suffix in order to make samples unique.
#' @param run.id.append When should run.ids be appended? If "all", all samples will have a run.id appeneded; if
#' "only.duplicated", only duplicated samples will have a run.id appended; if "none", run.ids will not be used to
#' make sample names unique.
#' @param duplicate.action The action to take to handle duplicates between marker genotype data frames. If \code{"merge"},
#' duplicates will be merged into a single consensus individual; if \code{"keep"}, duplicates will be retained in the final
#' output data frame; if \code{"keep.lowest.missing"}, the duplicate with the lowest missing data rate will be retained.
#'
#' @return A data.frame of merged genotype calls, where the first three columns are "marker," "chrom,", and "position."
#'
#' @importFrom dplyr coalesce
#'
#' @export
#'
merge_geno <- function(geno.list, run.ids = NULL, sep = ".", run.id.append = c("all", "only.duplicated", "none"),
                       duplicate.action = c("merge", "keep", "keep.lowest.missing")) {

  stopifnot(is.list(geno.list))
  stopifnot(is.character(sep))
  run.id.append <- match.arg(run.id.append)

  for (geno in geno.list) {
    stopifnot(is.data.frame(geno))
    cols.match <- names(geno)[1:3] == c("marker", "chrom", "position")
    if (any(!cols.match)) {
      stop("The first 3 columns of 'geno' should be 'marker', 'chrom', and 'position'.")
    }
  }

  if (!is.null(run.ids)) {
    stopifnot(is.character(run.ids))
    stopifnot(length(run.ids) == length(geno.list))
  } else {
    run.ids <- replicate(length(geno.list), expr = NULL, simplify = FALSE)
  }
  names(geno.list) <- run.ids

  # Message to user
  geno.list.length <- sapply(geno.list, dim, simplify = FALSE)

  cat(sub(pattern = "X", replacement = length(geno.list), x = "\nNumber of marker genotype data frames in the input list: X."))
  cat("\nSummary of input data:")
  cat("\nRun ID\tNumber of SNPs\tNumber of individuals")

  for (j in seq_along(geno.list.length)) {
    msg <- "\nX\tY\tZ"
    msg <- sub(pattern = "X", replacement = names(geno.list.length)[j], x = msg)
    msg <- sub(pattern = "Y", replacement = geno.list.length[[j]][1], x = msg)
    msg <- sub(pattern = "Z", replacement = geno.list.length[[j]][2], x = msg)
    cat(msg)
  }


  complements = c("A" = "T", "C" = "G", "T" = "A", "G" = "C")

  # Get marker names
  marker_names_list <- lapply(X = geno.list, FUN = function(x) x$marker)
  # Get unique
  marker_names <- sort(unique(unlist(marker_names_list)))

  sample_names_list <- lapply(X = geno.list, FUN = function(x) colnames(x)[-1:-3])

  # Index of where the sample names will be in the larger GT matrix
  sample_idy_list <- list()
  for (i in seq_along(sample_names_list)) {
    if (i == 1) {
      sample_idy_list[[i]] <- seq_along(sample_names_list[[i]])
    } else {
      sample_idy_list[[i]] <- max(sample_idy_list[[i-1]]) + seq_along(sample_names_list[[i]])
    }
  }


  # Make samples unique
  if (run.id.append == "none") {
    sample_names_uniq <- make.unique(names = unlist(sample_names_list), sep = sep)
  } else if (run.id.append == "all") {
    sample_names_uniq <- mapply(sample_names_list, run.ids, FUN = function(x, y) cbind(x, y))
    sample_names_uniq <- do.call(rbind, sample_names_uniq)
    sample_names_uniq <- apply(X = sample_names_uniq, MARGIN = 1, FUN = paste0, collapse = sep)
  } else {
    dup_samples <- table(unlist(sample_names_list))
    dup_samples <- names(dup_samples[dup_samples > 1])
    sample_names_list1 <- mapply(sample_names_list, run.ids, FUN = function(x, y) {
      x[x %in% dup_samples] <- paste0(x[x %in% dup_samples], sep, y)
      x
    })
    sample_names_uniq <- unlist(sample_names_list1)
  }


  # Create an empty marker matrix
  gt <- matrix(as.numeric(NA), nrow = length(marker_names), ncol = length(sample_names_uniq),
               dimnames = list(marker_names, sample_names_uniq))

  # Iterate over each vcf
  for (i in seq_along(geno.list)) {
    geno <- geno.list[[i]]
    # convert to a matrix
    geno_mat_i <- rrblup2genomat(x = geno, transpose = FALSE)
    mar_names <- row.names(geno_mat_i)
    # Add it to the larger GT
    idy <- sample_idy_list[[i]]
    gt[mar_names, idy] <- geno_mat_i
  }

  # Get marker, chrom, and position information and merge it
  snp_info_merged <- lapply(geno.list, "[", 1:3)
  snp_info_merged <- do.call(rbind, snp_info_merged)
  snp_info_merged <- unique(snp_info_merged)
  snp_info_merged <- snp_info_merged[order(snp_info_merged$chrom, snp_info_merged$position), ]

  gt1 <- gt[snp_info_merged$marker, , drop = FALSE]
  # colnames(gt1) <- make.unique(colnames(gt1))

  # Merge the marker metadata with the genotype call matrix to create the final data frame
  geno_full <- as.data.frame(snp_info_merged)
  geno_full <- cbind(geno_full, gt1[geno_full$marker, , drop = FALSE])

  ## Deal with duplicates
  # if action = "keep", do nothing
  if (duplicate.action == "keep") {
    cat("\n\nNo duplicates have been removed from the data.")
    return(geno_full)

  } else {
    # First create a data frame of original sample IDs and consensus sample IDs
    ids <- colnames(gt)
    # Remove run ID
    ids1 <- ids
    for (runid in run.ids) {
      ids1 <- sub(pattern = paste0(".", runid), replacement = "", x = ids1)
    }
    alias_df <- data.frame(original_id = colnames(gt), consensus_ids = ids1)

    if (duplicate.action == "keep.lowest.missing") {
      geno_full1 <- drop_duplicates(geno = geno_full, alias = alias_df, keep = "lowest.missing")

    } else {
      ## Merge genotypes
      # First identify duplicate sets
      alias_df_dup <- alias_df[alias_df[[2]] %in% alias_df[[2]][duplicated(alias_df[[2]])], , drop = FALSE]

      # Subset geno_full for non-duplicates
      geno_full1 <- geno_full[, setdiff(names(geno_full), alias_df_dup[[1]])]

      # Split alias_df_dup
      alias_df_dup_split <- split(alias_df_dup, alias_df_dup[[2]])
      geno_dup_consensus <- NULL
      # Iterate over these duplicates
      for (j in seq_along(alias_df_dup_split)) {
        alias_df_dup_j <- alias_df_dup_split[[j]]
        alias_df_dup_j_geno <- geno_full[, alias_df_dup_j[[1]]]
        # Order columns from left to right in decreasing missing data order
        alias_df_dup_j_geno <- alias_df_dup_j_geno[, order(colMeans(is.na(alias_df_dup_j_geno)))]
        # Coalesce
        consensus_geno <- coalesce(!!!as.list(alias_df_dup_j_geno))
        consensus_geno <- as.matrix(consensus_geno)
        colnames(consensus_geno) <- alias_df_dup_j[[2]][1]

        # Bind
        geno_dup_consensus <- cbind(geno_dup_consensus, consensus_geno)
      }

      # Bind the consensus genotypes with the full genotypes
      geno_full1 <- cbind(geno_full1, geno_dup_consensus)

      cat(sub(pattern = "X", replacement = ncol(geno_full) - 3, x = "\n\nNumber of individuals in the original 'geno' object: X"))
      cat(sub(pattern = "X", replacement = ncol(geno_full) - ncol(geno_full1), x = "\nNumber of duplicates merged: X"))
      cat(sub(pattern = "X", replacement = ncol(geno_full1) - 3, x = "\nNumber of individuals in the output 'geno' object: X"))

    }

    # Return
    return(geno_full1)

  }

}
