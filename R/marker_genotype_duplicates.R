#' Find and remove highly similar samples based on markers
#'
#' @description
#' Compute the number of matching marker genotypes between pairs of samples. Return the
#' pairs of individuals with similarity above a threshold.
#'
#' @param geno A data.frame of genotype calls, where the first three columns are "marker," "chrom,", and "position."
#' Subsequent columns are genotype calls for individuals. This can also be the output of \code{compare_geno} in
#' case you want to set a different threshold without running the marker comparison.
#' @param thresh The threshold minimum proportion of matching marker genotypes to return pairs of individuals
#'
#' @returns
#'
#' \itemize{
#'    \item{compare_geno}: A data frame with pairs of individuals and their similarity (only above the provided threshold)
#'    \item{group_duplicates}: A numeric vector of group numbers, with element names equal to individuals
#'  }
#'
#' @examples
#' geno <- replicate(3000, rbinom(200, 2, 0.5))
#' dimnames(geno) <- list(paste0("gen", seq_len(nrow(geno))), paste0("marker", seq_len(ncol(geno))))
#'
#' out <- compare_geno(geno, thresh = 0.4)
#'
#' # Show potential duplicates
#' out
#'
#' # Refilter using a different threshold
#' out1 <- compare_geno(out, thresh = 0.41)
#' out1
#'
#' @export
#'
compare_geno <- function(geno, thresh = 0.95) {
  stopifnot(inherits(geno, c("data.frame", "geno.comparison")))
  stopifnot(thresh >= 0 & thresh <= 1)

  if (!inherits(geno, "geno.comparison")) {
    cols.match <- names(geno)[1:3] == c("marker", "chrom", "position")
    if (any(!cols.match)) {
      stop("The first 3 columns of 'geno' should be 'marker', 'chrom', and 'position'.")
    }
    geno <- rrblup2genomat(x = geno)
    # Compute the proportions
    data <- geno
    z <- compare_geno_cpp(data)
    id <- row.names(geno)
    dimnames(z) <- list(id, id)
    diag(z) <- as.numeric(NA)

  } else {
    z <- attr(geno, "pairwise.comparison")

  }

  above_thresh <- which(z > thresh, arr.ind = TRUE)
  dups <- data.frame(geno1 = row.names(z)[above_thresh[,1]], geno2 = colnames(z)[above_thresh[,2]],
                     pMatching = z[above_thresh], row.names = NULL)
  row.names(above_thresh) <- NULL

  # add the matching proportion and threshold to the dups output
  attr(dups, "pairwise.comparison") <- z
  attr(dups, "threshold") <- thresh
  attr(dups, "above.threshold") <- above_thresh
  class(dups) <- c(class(dups), "geno.comparison")

  return(dups)

}

#'
#' @param dups The output of \code{compare_geno}.
#'
#' @rdname compare_geno
#'
#' @export
#'
group_duplicates <- function(dups) {
  # Get unique ids
  ids <- unique(dups$geno1)
  grouped_ids <- setNames(rep(0, length(ids)), ids)

  group <- 1
  while (any(grouped_ids == 0)) {
    id_no_group <- names(which(grouped_ids == 0)[1])
    ids_in_group <- dups[dups$geno1 == id_no_group,2]
    ids_grouped <- union(id_no_group, ids_in_group)
    grouped_ids[match(x = ids_grouped, names(grouped_ids))] <- group

    group <- group + 1
  }

  return(grouped_ids)

}

#' Compare and drop potential duplicates in a marker genotype matrix
#'
#' @param geno A data.frame of genotype calls, where the first three columns are "marker," "chrom,", and "position."
#' Subsequent columns are genotype calls for individuals.
#' @param alias A data frame where the first column is the original individual name
#' and the second column is the new consensus individual name.
#' @param keep The preference for keeping a consensus individual. If \code{"assigned"}, the
#' entry in the first column of \code{alias} will be retained. If \code{"lowest.missing"}, the
#' entry among those in the second column of \code{alias} with the lowest missing data
#' proportion will be retained.
#'
#' @details
#' For \code{drop_duplicates}, the function will coalesce marker genotypes for samples
#' with the same consensus name. Individuals not found in \code{alias} will be returned
#' as is.
#'
#' @returns A marker genotype matrix where rows are individuals and columns are markers.
#'
#'
#' @export
#'
drop_duplicates <- function(geno, alias, keep = c("assigned", "lowest.missing")) {
  stopifnot(is.data.frame(geno))
  cols.match <- names(geno)[1:3] == c("marker", "chrom", "position")
  if (any(!cols.match)) {
    stop("The first 3 columns of 'geno' should be 'marker', 'chrom', and 'position'.")
  }
  stopifnot(is.data.frame(alias))
  keep <- match.arg(keep)

  geno_orig <- geno
  geno <- rrblup2genomat(x = geno)

  # First find the individuals that are NOT in the alias df
  not_idx <- which(! row.names(geno) %in% alias[[1]])
  # Put this aside
  geno_not_idx <- geno[not_idx, , drop = FALSE]

  # Split by consensus genotype
  alias_split <- split(alias, alias[[2]])

  # First find entries in alias that have a one-to-one consensus genotype (i.e. nothing to merge)
  idx_singletons <- which(sapply(alias_split, nrow) == 1)
  alias_singletons <- names(alias_split[idx_singletons])
  idx_singletons <- which(row.names(geno) %in% alias_singletons)
  # Set these aside as well
  geno_singletons <- geno[idx_singletons, , drop = FALSE]

  # For non-singletons, return the individual with the least amount of missing data
  idx_duplicates <- which(sapply(alias_split, nrow) > 1)
  alias_split_duplicates <- alias_split[idx_duplicates]
  duplicate_keep <- vector("character", length(alias_split_duplicates))

  for (i in seq_along(alias_split_duplicates)) {
    dat <- alias_split_duplicates[[i]]
    if (keep == "assigned") {
      dup_keep <- unique(dat[[2]])

    } else {
      idx <- match(dat[[1]], row.names(geno))
      geno_id_missingness <- rowMeans(is.na(geno[idx, , drop = FALSE]))
      dup_keep <- names(geno_id_missingness)[which.min(geno_id_missingness)]

    }

    duplicate_keep[i] <- dup_keep

  }

  # Find the rows with the consensus duplicates
  idx_dup_keep <- which(row.names(geno) %in% duplicate_keep)
  geno_keep <- geno[idx_dup_keep, , drop = FALSE]

  # Merge geno mats
  geno_new <- rbind(geno_keep, geno_singletons, geno_not_idx)
  # Sort
  geno_new <- geno_new[order(row.names(geno_new)), , drop = FALSE]

  # Print messages
  cat("\nNumber of individuals in the original 'geno' object:", nrow(geno))
  cat("\nNumber of duplicates removed:", nrow(geno) - nrow(geno_new))
  cat("\nNumber of individuals in the output 'geno' object:", nrow(geno_new))

  geno_output <- geno_orig[,c(1:3, which(colnames(geno_orig) %in% row.names(geno_new)))]

  return(geno_output)

}





#' Reduce a marker matrix by duplicates
#'
#' @description
#' Takes a marker matrix and a data frame of geno_ids and aliases
#' and returns the geno_id or alias with the lowest missing marker data rate
#'
#'
#' @param geno A marker genotype matrix, where rows are markers and columns are individuals.
#' @param alias A data frame of individual and alias information, where the first column
#' is the desired genotype ID, and the second column contains alias information, if applicable
#' @param sep The delimiter separating multiple aliases in the second column of \code{alias}.
#'
#' @export
#'
find_alias_geno_missing <- function(geno, alias, sep = ", ") {

  # Error handling
  stopifnot(is.matrix(geno))
  stopifnot(is.data.frame(alias))
  # Rename the columns
  colnames(alias)[1:2] <- c("geno_id", "synonyms")
  stopifnot(is.character(sep))

  # Assign groups to alias
  alias1 <- alias
  alias1$group <- as.numeric(as.factor(alias$geno_id))

  # Split by group
  alias1_split <- split(alias1, alias1$group)

  # Calculate missing data rate for each genotype
  geno_ind_missing <- colMeans(is.na(geno))

  # Iterate over alias_split
  alias_missing_list <- list()

  for (i in seq_along(alias1_split)) {
    alias_df <- alias1_split[[i]]

    # Get the geno_id and aliases
    ind1 <- alias_df$geno_id
    ind2 <- strsplit(x = alias_df$synonyms, split = sep)[[1]]
    inds <- as.character(na.omit(c(ind1, ind2)))

    res <- data.frame(group = unique(alias_df$group), geno_id = ind1, alias_min_missing = inds, p_missing = geno_ind_missing[inds], row.names = NULL)
    res$which_min <- FALSE
    res$which_min[which.min(res$p_missing)] <- TRUE

    alias_missing_list[[i]] <- res

  }

  # Keep the list
  # Create a data.frame with each geno_id and the alias with the lowest amount of missing data
  alias_missing_df <- do.call("rbind", alias_missing_list)
  alias_missing_df <- alias_missing_df[alias_missing_df$which_min, ]
  alias_missing_df <- alias_missing_df[,c("geno_id", "alias_min_missing", "p_missing")]

  # Add the synonyms back and edit
  alias_updated <- merge(alias_missing_df, alias)
  idx <- which(alias_updated$geno_id != alias_updated$alias_min_missing)
  for (i in idx) {
    geno_id_use <- alias_updated$alias_min_missing[i]
    geno_id_old <- alias_updated$geno_id[i]
    synonyms_use <- union(setdiff(strsplit(x = alias_updated$synonyms[i], split = ", ")[[1]], geno_id_use), geno_id_old)
    alias_updated$geno_id[i] <- geno_id_use
    alias_updated$synonyms[i] <- paste0(synonyms_use, collapse = sep)

  }

  alias_updated <- alias_updated[,c("geno_id","synonyms")]

  output <- list(alias_updated = alias_updated, alias_missing = alias_missing_df, alias_missing_list = alias_missing_list)

  return(output)

}














