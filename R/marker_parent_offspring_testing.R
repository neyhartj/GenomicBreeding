#' Conduct marker quality control using parent-offspring tests
#' 
#' @description 
#' Uses heuristics and expected genotype probabilities to determine aberrant markers and (optionally) correct them.
#' 
#' @param x A \code{marker.matrix} object.
#' @param family.list A list of data frames that each include the columns "sample" (sample name), "pedigree", "parent1", "parent2", and "f_gen" (inbreeding generation). If the elements also contain the columns "parent1_name" and "parent2_name", missing parent names will be converted to these names.
#' @param min.fam.size The minimum family size necessary to correct offspring/parent genotype states when offspring are polymorphic.
#' @param correct Logical: should marker genotypes be corrected? Correction only happens when possible.
#' @param mismatch.tol
#' 
#' @import progress
#' @importFrom dplyr bind_rows
#' 
#' @export
#' 
parent_offspring_test_mf <- function(x, family.list, min.fam.size = 5, correct = TRUE, mismatch.tol = 3) {
  
  # Error checking
  stopifnot(is.list(family.list))
  stopifnot(inherits(x, "marker.matrix"))
  stopifnot(is.numeric(min.fam.size))
  stopifnot(is.logical(correct))
  stopifnot(is.numeric(mismatch.tol))
  
  # Make sure family list contains data.frames
  if (any(!sapply(X = family.list, is.data.frame))) stop ("Elements of family.list must be data.frames.")
  # Make sure these data.frames contains the right columns
  if (any(!sapply(X = family.list, function(df) all( c("sample", "pedigree", "parent1", "parent2", "f_gen") %in% names(df) ) ))) {
    stop("Elements of family.list must contain the following columns: 'sample', 'parent1', 'parent2', 'f_gen'.")
  }
  
  marker.matrix <- x
  
  # Create a list for output
  family_list_out <- list()
  
  # Initialize a progress bar
  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)",
                         total = length(family.list), complete = "=", incomplete = " ", show_after = 3)
  
  # Iterate over each element in family.list
  for (i in seq_along(family.list)) {
    
    # Update the progress bar
    pb$tick()
    
    family <- family.list[[i]]
    pedigree <- as.character(unique(family$pedigree))
    
    # Subset samples from the marker matrix
    sample_markers <- marker.matrix[family$sample,,drop = FALSE]
    parents <- unlist(unique(family[c("parent1", "parent2")]))
    parents_names <- unlist(unique(family[names(family) %in% c("parent1_name", "parent2_name")]))
    parents_names <- if (is.null(parents_names)) names(parents) else parents_names
    par_isNA <- is.na(parents)
    
    # Subset parents
    parent_markers <- marker.matrix[na.omit(parents),,drop = FALSE]
    # Add NA vector for missing parents
    missing_parent_replace <- parents_names[par_isNA]
    missing_parent_marker_matrix <- matrix(NA, nrow = length(missing_parent_replace), ncol = ncol(parent_markers), 
                                           dimnames = list(missing_parent_replace, colnames(parent_markers)))
    
    # Combine
    parent_markers <- rbind(parent_markers, missing_parent_marker_matrix)
    
    # Run the test
    po_test_out <- parent_offspring_test(parents = parent_markers, offspring = sample_markers, min.fam.size = min.fam.size,
                                         f.gen = min(family$f_gen), correct = correct, mismatch.tol = mismatch.tol)
    
    # Calculate pre/post missingness
    missingness_df <- data.frame(
      population = rep(c("parent", "offspring"), 2), 
      test = rep(c("pre", "post"), each = 2),
      missingness = c(mean(is.na(parent_markers)), mean(is.na(sample_markers)), 
                      mean(is.na(po_test_out$corrected_genos$parents)), mean(is.na(po_test_out$corrected_genos$offspring)))
    )
    
    # Return a list of the results
    family_list_out[[pedigree]] <- c(po_test_out, list(missingness = missingness_df))
    
  } # Close the loop
  
  # Pull out the corrected genotypes
  family_corrected_genos <- lapply(X = family_list_out, "[[", "corrected_genos")
  parent_corrected_genos <- do.call("rbind", lapply(X = family_corrected_genos, "[[", "parents"))
  
  # Combine parents and offspring
  family_corrected_genos <- rbind2(
    parent_corrected_genos,
    do.call("rbind", lapply(X = family_corrected_genos, "[[", "offspring"))
  )
  # Convert to marker matrix; resolve duplicates
  class(family_corrected_genos) <- union("marker.matrix", class(family_corrected_genos))
  family_corrected_genos <- resolve_duplicate_geno(x = family_corrected_genos)
  
  # Find identical samples within family_corrected_genos
  indentical_samples1 <- duplicated(x = family_corrected_genos, MARGIN = 1)
  family_corrected_genos_identical <- family_corrected_genos[indentical_samples1,,drop = FALSE]
  family_corrected_genos1 <- family_corrected_genos[!indentical_samples1,,drop = FALSE]
  
  # Find those entries in marker.matrix that were not in family_corrected_genos
  marker.matrix1 <- marker.matrix[! row.names(marker.matrix) %in% row.names(family_corrected_genos1),,drop = FALSE]
  # Test if any parents in family_corrected_genos are duplicates of remaining entries in marker.matrix 
  # (offspring are not likely to be duplicated)
  
  parent_genos_test <- family_corrected_genos1[row.names(family_corrected_genos1) %in% unique(row.names(parent_corrected_genos)), ,drop = FALSE]
  identical_samples2 <- duplicated(x = rbind(parent_genos_test, marker.matrix1), MARGIN = 1)
  identical_samples2 <- names(identical_samples2)[identical_samples2]
  
  # Get the names of duplicates
  family_corrected_genos_identical2 <- rbind2(
    parent_genos_test[row.names(parent_genos_test) %in% identical_samples2,, drop = FALSE],
    marker.matrix1[row.names(marker.matrix1) %in% identical_samples2,, drop = FALSE],
    
  )
  
  # Full identical genotypes list
  all_identitical_genos <- rbind2(family_corrected_genos_identical, family_corrected_genos_identical2)
  
  # Combine the unique genotypes
  combined_unique_genos <- rbind2(family_corrected_genos1, marker.matrix1[! row.names(marker.matrix1) %in% identical_samples2,, drop = FALSE])
  
  ## Pull out the PO test statistics
  po_test_notes <- lapply(X = family_list_out, "[[", "notes")
  po_test_notes1 <- lapply(X = po_test_notes, as.data.frame, responseName = "count")
  po_test_notes_df <- mapply(po_test_notes1, names(po_test_notes1), FUN = function(.x, .y) cbind(.x, family = .y), SIMPLIFY = FALSE)
  po_test_notes_df <- do.call("rbind", po_test_notes_df)
  row.names(po_test_notes_df) <- NULL
  
  # Pull out missingness statistics
  po_test_missing <- lapply(X = family_list_out, "[[", "missingness")
  po_test_missing1 <- mapply(po_test_missing, names(po_test_missing), FUN = function(.x, .y) cbind(.x, family = .y), SIMPLIFY = FALSE)
  po_test_missing_df <- do.call("rbind", po_test_missing1)
  row.names(po_test_missing_df) <- NULL
  
  
  # Assign class
  class(combined_unique_genos) <- union("marker.matrix", class(combined_unique_genos))
  class(all_identitical_genos) <- union("marker.matrix", class(all_identitical_genos))
  
  
  # Return results
  list(
    marker.matrix = combined_unique_genos,
    duplicated.marker.matrix = all_identitical_genos,
    parent.offspring.test.notes = po_test_notes_df,
    parent.offspring.test.missingness = po_test_missing_df
  )
  
  
} # close the function





#' 
#' @describeIn parent_offspring_test_mf
#' Conduct a parent-offpring test on a single family.
#' 
#' @param parents A 2 x m matrix of parent marker genotypes.
#' @param offspring A n x m matrix of offspring marker genotypes.
#' @param f.gen The filial (inbreeding) generation of offspring in the family (i.e. \code{f.gen = 3} = F_3 offspring).
#' 
#' @export
#' 
parent_offspring_test <- function(parents, offspring, f.gen, correct = TRUE, min.fam.size = 5, mismatch.tol = 3) {
  
  ## Error checking
  stopifnot(is.matrix(parents))
  stopifnot(is.matrix(offspring))
  stopifnot(nrow(parents) == 2)
  
  # check contents of parent and offspring
  if (!all(parents %in% c(-1, 0, 1, NA))) stop("Parent genotypes must take the form {-1, 0, 1}")
  if (!all(offspring %in% c(-1, 0, 1, NA))) stop("Parent genotypes must take the form {-1, 0, 1}")
  
  stopifnot(ncol(parents) == ncol(offspring))
  stopifnot(is.logical(correct))
  
  stopifnot(is.numeric(mismatch.tol))
  stopifnot(is.numeric(min.fam.size))
  
  ## Create a data.frame of notes for each marker
  # Notes will include the following options:
  # - All offspring NA
  # - Parents polymorphic / offspring consistent
  # - Parents polymorphic / offspring inconsistent
  # - Parents monomorphic / offspring monomorphic and match
  # - Parents monomorphic / offspring monomorphic and do not match
  # - Parents monomorphic / offspring polymorphic
  # - All or some Parents NA / offspring monomorphic
  # - Some parents NA / offspring polymorphic
  # - All parents NA / offspring polymorphic
  
  # Number of markers
  nM <- ncol(parents)
  # Number of offspring
  nOff <- nrow(offspring)
  
  ## Calculate genotype probabilities given f.gen and parent genotypes
  geno_prob_use <- geno_prob(f.gen)
  
  ## Measure NA proportion in parents and offspring
  parent_NA <- colMeans(is.na(parents))
  offspring_NA <- colMeans(is.na(offspring))
  
  # Check poly/monomorphic for whole set
  parent_mono <- abs(colMeans(parents)) == 1
  offspring_mono <- abs(colMeans(offspring, na.rm = TRUE)) == 1
  
  ## Vector of genotype contingencies
  geno_conting <- setNames(seq(-1, 1), seq(-1, 1))
  
  ## Empty list to store outcome
  note_vector <- character(nM)
  
  ## Loop over markers
  for (j in seq_len(nM)) {
    
    # Pull out marker genotypes
    parent_j <- parents[,j]
    offspring_j <- offspring[,j]
    
    ## First check NA level of offspring
    # If all NA, note and continue
    if (offspring_NA[j] == 1) {
      
      # Check parents
      note_vector[j] <- ifelse(parent_NA[j] > 0, "Some or all parents NA/offspring NA",
                               "Parents observed/offspring NA")
      
    } else {
      # Check if parents are NA
      # If any NA, proceed
      if (parent_NA[j] > 0) {
        
        # Check if offspring are monomorphic
        if (offspring_mono[j]) {
          # Record this
          note_vector[j] <- "All or some Parents NA/offspring monomorphic"
          
          # Because we already checked for offspring NA above, if offspring are not
          # monomorphic, they must be polymorphic
          #
          # Check number of NA in parents; note this
        } else {
          note_vector[j] <- ifelse(parent_NA[j] == 0.5, "Some parents NA/offspring polymorphic",
                                   "All parents NA/offspring polymorphic")
          
        }
        
        # Now check if parents are monomorphic/polymorphic
      } else if (parent_mono[j]) {
        
        # Check if offspring are monomorphic
        if (offspring_mono[j]) {
          
          # Check if all non-missing offspring match the parents; note this
          note_vector[j] <- ifelse(all(na.omit(offspring_j) == unique(parent_j)),
                                   "Parents monomorphic/offspring monomorphic and match",
                                   "Parents monomorphic/offspring monomorphic and do not match")
          
          
          # else offspring are polymorphic
        } else {
          # Find the number of offspring that do not match the parents
          n_mismatch <- sum(na.omit(offspring_j) != unique(parent_j))
          
          # Report this
          note_vector[j] <- paste0("Parents monomorphic/offspring polymorphic mismatches: ", n_mismatch)
          
        }
        
        # Else parents polymorphic
      } else {
        
        # Calculate expected genotype propabilities of 0, 1, 2
        pGeno <- geno_prob_use[paste0(sort(parent_j), collapse = "|"),,drop = TRUE]
        
        ## Custom X2_test
        obs <- sapply(X = geno_conting, FUN = function(type) sum(offspring_j == type, na.rm = TRUE))
        exp <- pGeno * sum(obs)
        stat <- sum(((obs - exp)^2) / exp)
        p_value <- pchisq(q = stat, df = length(obs) - 1, lower.tail = FALSE)
        
        # Look at p-value, note
        note_vector[j] <- ifelse(p_value >= 0.05, "Parents polymorphic/offspring consistent",
                                 "Parents polymorphic/offspring inconsistent")
        
      }
      
    }
    
    # Close loop
  }
  
  ## Contingency table of the note vector; store vector of marker indices for each outcome
  note_table <- table(note_vector)
  
  ## Split marker indices by note type
  note_index_split <- split(seq_along(note_vector), note_vector)
  
  
  ## Corrections ##
  
  # Ony a few contigencies can be corrected absent phased genotypic data:
  # 1. All or some parents NA / offspring monomorphic: correct parents
  # 2. Parents monomorphic and offspring polymorphic with mismatch:
  # 2a. if below some mismatch threshold, hide the offspring genotype
  # 2b. if above some mismatch threshold, make parents NA
  # 3. Parents polymorphic/offspring inconsistent: choose parent genotype that maximizes likelihood
  # 4. Some parents NA/offspring polymorphic: correct missing parent
  
  if (correct) {
    
    # Create correction matrices
    parents_corrected <- parents
    offspring_corrected <- offspring
    
    ## Address contingencies one-at-a-time
    ##
    ## 1. All or some parents NA / offspring monomorphic: correct parents and any missing offspring
    index_to_correct <- note_index_split$`All or some Parents NA/offspring monomorphic`
    # subset
    parents_to_correct <- parents[,index_to_correct, drop = FALSE]
    offspring_to_correct <- offspring[,index_to_correct, drop = FALSE]
    
    # Determine unique offspring genotypes for each marker
    offspring_unique_geno <- apply(X = offspring_to_correct, MARGIN = 2, FUN = function(x) unique(na.omit(x)))
    
    # Skip if empty
    if (length(offspring_unique_geno) > 0) {
      
      ## Apply to parents
      parents_to_correct <- matrix(data = replicate(2, offspring_unique_geno), nrow = nrow(parents),
                                   ncol = ncol(offspring_to_correct), dimnames = list(row.names(parents), colnames(offspring_to_correct)),
                                   byrow = TRUE)
      ## Apply to offspring
      offspring_to_correct <- matrix(data = replicate(nrow(offspring_to_correct), offspring_unique_geno), nrow = nrow(offspring_to_correct),
                                     ncol = ncol(offspring_to_correct), dimnames = dimnames(offspring_to_correct), byrow = TRUE)
      
      ## Replace
      parents_corrected[,index_to_correct] <- parents_to_correct
      offspring_corrected[,index_to_correct] <- offspring_to_correct
      
    }
    
    
    # 2. Parents monomorphic and offspring polymorphic with mismatch:
    ## Find the mismatch vectors with mismatch below the threshold
    # index of list names with mismatch
    
    # 2a - correct offspring
    mismatch_index_2a <- which(as.numeric(str_extract(names(note_index_split), "[0-9]{1,}")) <= mismatch.tol)
    index_to_correct <- unlist(note_index_split[mismatch_index_2a])
    
    # subset
    parents_to_correct <- parents[,index_to_correct, drop = FALSE]
    offspring_to_correct <- offspring[,index_to_correct, drop = FALSE]
    
    ## Loop over markers
    for (j in seq_len(ncol(parents_to_correct))) {
      # What is the intended genotype state
      state_to_use <- unique(parents_to_correct[,j])
      # Find the mismatches (and NAs) in the offspring - replace
      offspring_to_correct[which(offspring_to_correct[,j] != state_to_use | is.na(offspring_to_correct[,j])), j] <-
        state_to_use
    }
    
    # Replace
    parents_corrected[,index_to_correct] <- parents_to_correct
    offspring_corrected[,index_to_correct] <- offspring_to_correct
    
    # 2b. Correct parents
    mismatch_index_2b <- which(as.numeric(str_extract(names(note_index_split), "[0-9]{1,}")) > mismatch.tol)
    index_to_correct <- unlist(note_index_split[mismatch_index_2b])
    
    # subset
    parents_to_correct <- parents[,index_to_correct,drop = FALSE]
    # Set all as NA
    parents_to_correct[,] <- NA
    
    # Replace
    parents_corrected[,index_to_correct] <- parents_to_correct
    
    
    ## 3. Parents polymorphic/offspring inconsistent: choose parent genotype that maximizes likelihood
    ## 4. Some parents NA/offspring polymorphic: correct missing parent
    index_to_correct <- unlist(note_index_split[c("Parents polymorphic/offspring inconsistent",
                                                  "Some parents NA/offspring polymorphic")])
    # subset
    parents_to_correct <- parents[,index_to_correct, drop = FALSE]
    offspring_to_correct <- offspring[,index_to_correct, drop = FALSE]
    
    ## Calculate new genotype probabilities; ignore hets combinations
    geno_prob_use1 <- geno_prob_use[c("-1|1", "-1|-1", "1|1"),]
    
    ## Calculate the parent genotypes that maximize the likelihood based on offspring
    ## segregation
    ## 
    ## Only do this if enough offspring are present
    ## 
    if (nOff >= min.fam.size) {
      parent_genotype_max_likelihood <- apply(X = offspring_to_correct, MARGIN = 2, FUN = function(x) {
        # Custom table calculation
        off_table <- sapply(X = geno_conting, FUN = function(type) sum(x == type, na.rm = TRUE))
        apply(X = geno_prob_use1, MARGIN = 1, FUN = dmultinom, x = off_table,
              size = sum(off_table), log = TRUE)
      })
      
    } else {
      parent_genotype_max_likelihood <- numeric()
      
    }
    
    if (length(parent_genotype_max_likelihood) > 0) {
      
      # Find the genotype states that maximize the logLik
      max_log_lik <- row.names(geno_prob_use1)[apply(X = parent_genotype_max_likelihood, MARGIN = 2, FUN = which.max)]
      # Split the codes
      max_log_lik1 <- lapply(str_split(string = max_log_lik, pattern = "\\|"), as.numeric)
      
      ## Iterate over parent states and correct, if necessary
      for (j in seq_len(ncol(parents_to_correct))) {
        parent_state <- parents_to_correct[,j]
        
        # Compare to the max loglik genotype states
        if (all(max_log_lik1[[j]] == sort(parent_state))) {
          # If they match, continue
          next
          
        } else {
          # 1. Are the replacements identical? If so, just replace all parents
          if (n_distinct(max_log_lik1[[j]]) == 1) {
            parents_to_correct[,j] <- max_log_lik1[[j]]
            
            # 2. Do both parents need to be replaced? If so, set both as NA (cannot differentiate)
          } else if (!any(parent_state %in% max_log_lik1[[j]])) {
            parents_to_correct[,j] <- NA
            
            # Else find the parent to be replaced
          } else {
            # Find the mismatch in the parents
            missing_state <- setdiff(max_log_lik1[[j]], parent_state)
            missing_state <- ifelse(is_empty(missing_state), unique(max_log_lik1[[j]]), missing_state)
            
            # Which parent needs correction
            parent_j1 <- which(! parent_state %in% max_log_lik1[[j]])
            parents_to_correct[parent_j1, j] <- missing_state
            
          }
        }
        
        ## Next if parent genotypes are the same, impute on offspring
        if (n_distinct(parents_to_correct[,j]) == 1) {
          offspring_to_correct[,j] <- unique(parents_to_correct[,j])
          
        }
        
      }
      
      ## Replace
      parents_corrected[,index_to_correct] <- parents_to_correct
      offspring_corrected[,index_to_correct] <- offspring_to_correct
      
    } # Close the if statement
    
  } else {
    # Else no corrections
    parents_corrected <- parents
    offspring_corrected <- offspring
    
  }
  
  
  ## Output a list
  out <- list(corrected_genos = list(parents = parents_corrected, offspring = offspring_corrected),
              notes = note_table)
  
  return(out)
  
}




