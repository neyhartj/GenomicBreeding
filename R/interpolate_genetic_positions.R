#' Predict marker genetic position using linear interpolation
#'
#' @description
#' Uses linear interpolation and local recombination rate information to predict
#' the unknown genetic position of markers based on their known physical positions.
#'
#' @param map A data frame with three columns: marker name, chromosome, and physical
#' position. These are markers whose genetic position is unknown
#' @param cMMb A data frame with six columns: chromosome, left flank physical position,
#' right flank physical position, local recombination rate (in cM/Mb) in that flank,
#' left flank genetic position, and right flank genetic position.
#'
#'
#' @export
#'
interpolate_genetic_position <- function(map, cMMb) {
  # Error checking
  stopifnot(is.data.frame(map))
  stopifnot(is.data.frame(cMMb))

  # All chromosomes in map in cmmb?
  if (!all(unique(map[[2]]) %in% unique(cMMb[[1]]))) stop ("The chromosomes in 'map' are not all in 'cMMb'.")

  # Add a new column to map
  map1 <- map
  map1$gen_pos <- as.numeric(NA)

  # Iterate over markers in the map
  for (i in seq_len(nrow(map1))) {
    # Get the position
    pos_i <- map1$pos[i]
    chrom_i <- map1$chrom[i]
    # Find the interval this position belongs to
    interval_i <- which(cMMb[[1]] == chrom_i & cMMb[[2]] <= pos_i & cMMb[[3]] >= pos_i)
    # If empty, skip
    if (is_empty(interval_i)) next

    # Predict the genetic position of the marker from the start of the interval
    dist_i <- pos_i - cMMb[[2]][interval_i]
    dist_cm_i <- cMMb[[4]][interval_i] * (dist_i / 1e6)
    pred_gen_i <- cMMb[[5]][interval_i] + dist_cm_i

    # add this to the df
    map1$gen_pos[i] <- pred_gen_i

  }

  # Return the adjusted map
  return(map1)

}
