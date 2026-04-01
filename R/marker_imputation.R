#' Impute marker genotypes in families using flanking markers
#' 
#' @description 
#' Imputed marker genotypes for multiple families using flanking markers and recombination frequencies.
#' 
#' @param x A \code{marker.matrix} object.
#' @param family.list A list of data frames that each include the columns "sample" (sample name), "pedigree", "parent1", "parent2", and "f_gen" (inbreeding generation). If the elements also contain the columns "parent1_name" and "parent2_name", missing parent names will be converted to these names.
#' @param marker.map A data frame of marker map positions, where the first column is the marker name, the second column is the chromosome, and the third column is the genetic map position.
#' @param min.geno.prob The minimum marker genotype probability needed for imputation.
#' @param map.function The mapping function to convert cM distance to recombination frequency.
#' @param return.geno.prob Logical; should genotype probabilities also be returned?
#' 
#' @return 
#' A list containing i) the new marker matrix with imputed marker genotypes and ii) the missingness statistics for pre-  and post-imputation.
#' 
#' @import progress
#' 
#' @export
#' 
#' 
impute_geno_flanking_mf <- function(x, family.list, marker.map, min.geno.prob = 0.7, 
                                    map.function = c("haldane","kosambi","c-f","morgan"), return.geno.prob = FALSE) {
  
  # Error checking
  marker.matrix <- x
  stopifnot(is.list(family.list))
  stopifnot(inherits(marker.matrix, "marker.matrix"))
  stopifnot(is.data.frame(marker.map))
  stopifnot(is.logical(return.geno.prob))
  
  # Make sure all markers in marker.matrix are in the marker map
  if (!all(colnames(marker.matrix) %in% marker.map[[1]])) stop("Not all of the markers in marker.matrix are in the marker.map object.")
  
  # Sort the marker map
  marker.map1 <- marker.map[order(marker.map[[2]], marker.map[[3]]),]
  marker.map2 <- marker.map1[marker.map1[[1]] %in% colnames(marker.matrix),]
  
  # Copy the marker.matrix object to overlay imputed data; also copy to add genotype probabilities
  marker_geno_prob <- marker.matrix1 <- marker.matrix
  # Add NA probabilities
  marker_geno_prob[,] <- NA
  class(marker_geno_prob) <- setdiff(class(marker_geno_prob), "marker.matrix")
  
  
  # Make sure family list contains data.frames
  if (any(!sapply(X = family.list, is.data.frame))) stop ("Elements of family.list must be data.frames.")
  # Make sure these data.frames contains the right columns
  if (any(!sapply(X = family.list, function(df) all( c("sample", "pedigree", "parent1", "parent2", "f_gen") %in% names(df) ) ))) {
    stop("Elements of family.list must contain the following columns: 'sample', 'parent1', 'parent2', 'f_gen'.")
  }
  
  # If family.list has no names, set to the pedigree
  if (is.null(names(family.list))) {
    names(family.list) <- sapply(X = family.list, FUN = function(x) unique(x$pedigree))
  }
  
  # Create a data.frame for output
  imputation_stats_out <- data.frame(pedigree = names(family.list), missingness_pre = NA, missingness_post = NA, stringsAsFactors = FALSE)
  
  # Initialize a progress bar
  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)",
                         total = length(family.list), complete = "=", incomplete = " ", show_after = 3)
  
  # Iterate over each element in family.list
  for (i in seq_along(family.list)) {
    
    # Update the progress bar
    pb$tick()
    
    family <- family.list[[i]]
    f_gen = min(family$f_gen)
    
    # Offspring names
    offspring_names <- family$sample
    # Parent names
    parent_names <- unlist(unique(family[c("parent1", "parent2")]))
    
    # Are both parents in the marker matrix
    both_parents_in_marker_matrix <- all(parent_names %in% row.names(marker.matrix))
    # If not, skip this family
    if (!both_parents_in_marker_matrix) next
    
    # Subset marker genotypes for each set
    parents <- marker.matrix[parent_names, marker.map2[[1]], drop = FALSE]
    offspring <- marker.matrix[offspring_names, marker.map2[[1]], drop = FALSE]
    
    # Run the imputation function
    offspring_imputed_out <- impute_geno_flanking(parents = parents, offspring = offspring, marker.map = marker.map2, 
                                                  f.gen = f_gen, min.geno.prob = min.geno.prob, map.function = map.function,
                                                  return.geno.prob = return.geno.prob)
    
    # Replace the genotypes for these offspring in the marker.matrix object
    marker.matrix1[row.names(offspring_imputed_out), ] <- offspring_imputed_out
    
    # Add marker genotype probabilities
    if (return.geno.prob) {
      marker_geno_prob[row.names(offspring_imputed_out), ] <- attr(offspring_imputed_out, "geno.prob")
    }
    
    # Add the stats to imputation_stats_out
    imputation_stats_out[i,c("missingness_pre", "missingness_post")] <- c(mean(is.na(offspring)), mean(is.na(offspring_imputed_out)))
    
  } # Close the family loop
  
  # Reclass the marker matrix
  class(marker.matrix1) <- class(marker.matrix)

  # Return the marker matrix and the missingness stats
  if (return.geno.prob) {
    out <- list(
      marker.matrix.imputed = marker.matrix1,
      marker.geno.prob = marker_geno_prob,
      missingness.stats = imputation_stats_out
    )
    
  } else {
    out <- list(
      marker.matrix.imputed = marker.matrix1,
      missingness.stats = imputation_stats_out
    )
    
  }
  
  return(out)
    
} # Close the function




#' 
#' @describeIn impute_geno_flanking_mf
#' Impute marker genotypes in a single family using flanking markers
#' 
#' @param parents A 2 x m matrix of parent marker genotypes.
#' @param offspring A n x m matrix of offspring marker genotypes.
#' @param f.gen The filial (inbreeding) generation of offspring in the family (i.e. f.gen = 3 = F_3 offspring).
#' 
#' @return
#' A matrix of imputed offspring marker genotypes.
#' 
#' @importFrom qtl table2map sim.cross jittermap calc.genoprob
#' 
#' @export 
#' 
#' 
impute_geno_flanking <- function(parents, offspring, marker.map, f.gen, min.geno.prob = 0.7, 
                                 map.function = c("haldane","kosambi","c-f","morgan"), return.geno.prob = FALSE) {
  
  # Error checking
  stopifnot(is.matrix(parents))
  stopifnot(is.matrix(offspring))
  stopifnot(nrow(parents) == 2)
  
  # check contents of parent and offspring
  if (!all(parents %in% c(-1, 0, 1, NA))) stop("Parent genotypes must take the form {-1, 0, 1}")
  if (!all(offspring %in% c(-1, 0, 1, NA))) stop("Parent genotypes must take the form {-1, 0, 1}")
  
  stopifnot(ncol(parents) == ncol(offspring))
  stopifnot(identical(colnames(parents), colnames(offspring)))
  
  stopifnot(is.numeric(f.gen))
  stopifnot(is.numeric(min.geno.prob))
  
  # Convert marker.map to a qtl map object
  marker.map1 <- as.data.frame(marker.map[-1])
  row.names(marker.map1) <- marker.map[[1]]
  marker_map_use <- table2map(tab = marker.map1)
  
  # Index of all markers
  nM <- ncol(parents)
  marker_index <- seq_len(nM)
  
  # Rename
  parent_states <- parents
  offspring1 <- offspring
  
  # First find any markers that are monomorphic in the parents and impute on the offspring
  monomorphic_par_markers <- abs(colSums(parent_states, na.rm = TRUE)) == 2
  dimorphic_par_markers <- abs(colSums(parent_states, na.rm = TRUE)) == 0 & colSums(is.na(parent_states)) == 0
  # Get the single marker genotype of these monomorphic markers
  monomorphic_par_markers_geno <- parent_states[1,monomorphic_par_markers]
  
  # Impute on the offspring
  offspring1[,monomorphic_par_markers] <- matrix(monomorphic_par_markers_geno, nrow = nrow(offspring1), 
                                                 ncol = sum(monomorphic_par_markers), byrow = TRUE)
  
  # Next find the index of markers where there are any missing offspring
  any_missing_off_markers <- which(colSums(is.na(offspring1)) > 0)
  
  # Assign parent states to the markers
  parent_states[1,!is.na(parent_states[1,])] <- 1
  parent_states[2,!is.na(parent_states[2,])] <- -1
  
  # In the offspring, only assign parent states that are dimorphic in the parents
  offspring_states <- offspring1
  for (j in seq_len(ncol(offspring_states))) {
    # Is the marker dimorphic?
    if (dimorphic_par_markers[j]) {
      # Assign parent states
      par_j <- parents[,j]
      off_j <- offspring1[,j]
      off_states_j <- ifelse(off_j == par_j[1], 1, -1)
      off_states_j[off_j == 0] <- 0
      
    } else {
      # Otherwise the state is unknown (NA)
      off_states_j <- NA
      
    }
    
    # Replace off state
    offspring_states[,j] <- off_states_j
    
  }
  
  # Next use R/qtl to calculate the genotype probabilities
  
  
  ## Create a cross object
  # First adjust the map
  marker_map_use1 <- jittermap(marker_map_use, 1e-10)
  nInd <- nrow(offspring)
  
  # Simulate a cross
  cross_new <- sim.cross(map = marker_map_use1, n.ind = nInd, type = "bcsft", 
                         cross.scheme = c(s = 0, t = f.gen))
  
  ## Overlay the genotypes into the cross data
  for (chr in seq_along(marker_map_use1)) {
    cross_new$geno[[chr]]$data <- offspring_states[,names(marker_map_use1[[chr]]), drop = FALSE] + 2
  }
  
  ## Calculate genotype probabilities
  cross1_genoprob <- calc.genoprob(cross = cross_new, step = 0, off.end = 0, 
                                   map.function = map.function, stepwidth = "fixed")
  
  # Extract the genotype probabilities
  cross1_genoprob1 <- lapply(X = cross1_genoprob$geno, FUN = "[[", "prob")
  # Get the genotype with the highest probability
  cross1_most_prob_geno <- lapply(X = cross1_genoprob1, FUN = function(ar) apply(X = ar, MARGIN = c(1,2), FUN = which.max) )
  cross1_most_prob_geno <- do.call("cbind", cross1_most_prob_geno)
    
  # Get the probability of the most probable genotype
  cross1_prop_most_prob <- lapply(X = cross1_genoprob1, FUN = function(ar) apply(X = ar, MARGIN = c(1,2), FUN = max) ) 
  cross1_prop_most_prob <- do.call("cbind", cross1_prop_most_prob)
  
  # Any parent state genotypes with probabilities < min.geno.prob set to NA
  offspring_states1 <- cross1_most_prob_geno
  offspring_states1[cross1_prop_most_prob < min.geno.prob] <- NA
  
  
  
  
  # Run through offspring and overlay parent genotypes (not states) on the offspring given the most probable genotype
  # 
  
  # Create a matrix of parent1 / het / parent2
  parents1 <- rbind(parents[1,,drop = FALSE], het = 0, parents[2,,drop = FALSE])
  
  # Overlay the parent genotypes
  offspring_imputed_genotypes <- offspring1
  # Iterate
  for (j in seq_len(ncol(offspring_imputed_genotypes))) {
    # Which offspring are NA for this SNP
    na_j <- is.na(offspring_imputed_genotypes[,j])
    offspring_imputed_genotypes[na_j,j] <- parents1[offspring_states1[na_j,j],j]
  }
  
  # Add the genotype probabilities, if called
  if (return.geno.prob) {
    attr(offspring_imputed_genotypes, "geno.prob") <- cross1_prop_most_prob
    
  }
  
  # Return the imputed offspring genotypes
  return(offspring_imputed_genotypes)
  
}




#' 
#' @describeIn impute_geno_flanking_mf
#' Impute marker genotypes using an expectation-maximization algorithm or the marker mean; implements this one chromosome at-a-time.
#' 
#' @param method The imputation method; may be "EM" for expectation-maximization or "mean" for the marker mean.
#' @param verbose Logical; should status messages be printed?
#' 
#' @importFrom rrBLUP A.mat
#' 
#' @export
#' 
#' 
impute_geno_em <- function(x, method = c("EM", "mean"), verbose = TRUE) {
  
  # Error checking
  stopifnot(inherits(x, "marker.matrix"))
  method <- match.arg(method)
  stopifnot(is.logical(verbose))
  
  # Get the snp metadata
  snp_meta <- attr(x, "snp.metadata")
  snp_cols <- attr(snp_meta, "cols")
  
  # Get a list of markers per chromosome
  snp_meta_split <- split(x = snp_meta, f = snp_meta[[snp_cols[["chrom.col"]]]])
  snps_per_chrom <- lapply(X = snp_meta_split, FUN = "[[", snp_cols[["name.col"]])
  
  # Create a list to store the imputed marker
  imputed_out_list <- vector("list", length = length(snps_per_chrom))
  names(imputed_out_list) <- names(snps_per_chrom)
  
  # Iterate over the chromosomes
  for (i in seq_along(snps_per_chrom)) {
    # Names of snps in this chromosome
    snps_chrom_i <- snps_per_chrom[[i]]
    # Subset the marker matrix
    x_i <- x[, snps_chrom_i, drop = FALSE]
    
    # Implement the em algorithm
    temp <- capture.output(
      impute_out <- A.mat(X = x_i, min.MAF = 0, max.missing = 1, impute.method = method, return.imputed = TRUE)
    )
    
    # Print a message if verbose
    if (verbose) {
      quotes <- gregexpr(pattern = '"', text = temp)[[1]]
      temp_sub <- substr(x = temp, start = quotes[1] + 1, stop = quotes[2] - 1)
      msg <- paste0("\nImputation of chromosome ", i, " complete.\n", temp_sub, "\n")
      cat(msg)
      
    }
    
    # Add the imputed marker matrix to the list
    imputed_out_list[[i]] <- impute_out$imputed
    
  }
  
  # Combine the imputed genotypes
  imputed_out_combined <- do.call("cbind", imputed_out_list)
  
  # Give it the same attributes
  attributes(imputed_out_combined) <- attributes(x)

  # Return
  return(imputed_out_combined)
  
}





