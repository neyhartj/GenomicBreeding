#' Read or write marker genotype data
#' 
#' @description 
#' Reads marker genotype data to a \code{marker.geno} object, or writes marker genotype data to a file.
#' 
#' @param x A \code{data.frame} or \code{marker.matrix} to be converted into a \code{marker.geno} object.
#' 
#' @importFrom dplyr as_tibble
#' 
#' @export 
#' 
as_marker_geno <- function(x, ...) {
  UseMethod("as_marker_geno", x)
}

#' 
#' @rdname as_marker_geno
#' 
#' @param chrom.col The name of the column with chromosome information. (\emph{Required})
#' @param pos.col The name of the column with physical chromosome positions. (\emph{Required})
#' @param cMpos.col The name of the column with centi-Morgan positions.
#' @param allele.col The name of the column with allele information.
#' 
#' @export
#'  
as_marker_geno.data.frame <- function(x, chrom.col, pos.col, cmPos.col = NULL, allele.col = NULL) {
  
  # Error checking
  stopifnot(is.data.frame(x))
  stopifnot(is.character(chrom.col))
  stopifnot(is.character(pos.col))
  if (!is.null(cmPos.col)) stopifnot(is.character(cmPos.col)) else cmPos.col <- NULL
  if (!is.null(allele.col)) stopifnot(is.character(allele.col)) else allele.col <- NULL
  
  # Rename x
  data_in <- x
  
  # Get the name of the marker name column
  name.col <- names(data_in)[1]
  # Vector of columns to avoid
  meta_cols <- c(name.col, chrom.col, pos.col, cmPos.col, allele.col)
  meta_cols <- meta_cols[!is.null(meta_cols)]
  # Subset the data without metadata
  data_in1 <- data_in[! names(data_in) %in% meta_cols]
  
  # Detect the genotype format (AA/BB; ATCG; 012; -101)
  geno_formats <- sapply(X = data_in1, FUN = detect_genotype_format)
  geno_formats_prop <- prop.table(x = table(geno_formats))
  
  # If the genotype formats are < 95% for a single format, error out
  if (max(geno_formats_prop) < 0.75 | names(which.max(geno_formats_prop)) %in% "unclear") {
    stop ("The genotype formats in the file are ambiguous; check them for inconsistent formats.")
    
  } else {
    geno_format <- names(which.max(table(geno_formats)))
    
  }
  
  ## Convert the data to -101 ##
  data_converted <- as.matrix(data_in1)
  
  # First error if the format is "ATCG" and no allele data is present
  if (geno_format == "ATCG" & is.null(allele.col)) 
    stop("The genotype format cannot be nucleotide (i.e. 'ATCG') with a missing allele column.")
  
  # Convert
  if (geno_format == "ATCG") {
    stop("Nucleotide format is not yet supported.")
    
  } else if (geno_format == "AA/BB") {
    data_converted[data_converted == "AA"] <- 1
    data_converted[data_converted == "AB"] <- 0
    data_converted[data_converted == "BB"] <- -1
    
    data_converted1 <- apply(X = data_converted, MARGIN = 2, FUN = parse_number, 
                             na = setdiff(unique(as.vector(data_converted)), c("-1", "0", "1")))
    
    # Convert back to data.frame
    data_in2 <- as.data.frame(data_converted1)
    
  } else if (geno_format == "012") {
    # Parse to number
    if (is.character(data_converted)) {
      data_converted1 <- apply(X = data_converted, MARGIN = 2, FUN = parse_number, 
                               na = setdiff(unique(as.vector(data_converted)), c("0", "1", "2")))
      
    } else {
      data_converted1 <- data_converted
      
    }
    
    # Substract 1
    data_converted1 <- data_converted1 - 1
    
  } else if (geno_format == "-101") {
    # Parse to number
    if (is.character(data_converted)) {
      data_converted1 <- apply(X = data_converted, MARGIN = 2, FUN = parse_number, 
                               na = setdiff(unique(as.vector(data_converted)), c("0", "1", "2")))
      
    } else {
      data_converted1 <- data_converted
      
    }
    
  }
  
  # Convert back to data.frame
  data_converted2 <- as.data.frame(data_converted1)
  # Merge with the metadata
  data_in2 <- cbind(data_in[names(data_in) %in% meta_cols], data_converted2)
  
  # Vector of metadata columns
  meta_cols <- c("name.col" = name.col, "chrom.col" = chrom.col, "pos.col" = pos.col, 
                 "cmPos.col" = cmPos.col, "allele.col" = allele.col)
  
  ## Edit names
  meta_cols_pos <- which(names(data_in2) %in% meta_cols)
  # Capitalize entry names aside from the meta cols
  names(data_in2)[-meta_cols_pos] <- toupper(names(data_in2)[-meta_cols_pos])
  # Remove spaces
  names(data_in2)[-meta_cols_pos] <- gsub(pattern = " ", replacement = "-", x = names(data_in2)[-meta_cols_pos])
  # Make unique
  names(data_in2) <- make.unique(names = names(data_in2), sep = ";")
  
  # Identify the duplicates; store this information as an attribute in the df
  dups <- unique(gsub(pattern = ";[0-9]{1,}$", replacement = "", x = grep(pattern = ";[0-9]{1,}$", x = names(data_in2), value = TRUE)))
  dups <- if (length(dups) == 0) NULL else dups
  
  # Convert to tibble
  data_in3 <- as_tibble(data_in2)
  
  # Add other attributes
  class(data_in3) <- union("marker.geno", class(data_in3))
  attr(data_in3, "cols") <- meta_cols
  attr(data_in3, "duplicates") <- dups
  
  
  # Sort by chromosome, then position
  data_in4 <- data_in3[order(data_in3[[chrom.col]], data_in3[[pos.col]]),]
  # Return the object
  return(data_in4)
  
}

#'
#' @rdname as_marker_geno
#' 
#' @param chrom.vec A character vector of chromosome assignments; the length of this vector must equal the number of markers.
#' @param pos.vec A numeric vector of positions; the length of this vector must equal the number of markers.
#' @param cmPos.vec A numeric vector of genetic positions; the length of this vector must equal the number of markers.
#' @param allele.vec A character vector of marker alleles; the length of this vector must equal the number of markers.
#' 
#' @export
#'
as_marker_geno.marker.matrix <- function(x, chrom.vec, pos.vec, cmPos.vec = NULL, allele.vec = NULL) {
  
  # Error checking
  stopifnot(is.character(chrom.vec))
  stopifnot(is.numeric(pos.vec))
  stopifnot(all.equal(length(chrom.vec), length(pos.vec)))
  stopifnot(all.equal(ncol(x), length(chrom.vec)))
  
  if (!is.null(cmPos.vec)) {
    stopifnot(is.numeric(cmPos.vec))
    stopifnot(all.equal(length(cmPos.vec), length(chrom.vec)))
    
  }
  
  if (!is.null(allele.vec)) {
    stopifnot(is.character(allele.vec))
    stopifnot(all.equal(length(allele.vec), length(chrom.vec)))
    
  }
  
  ## Convert the marker matrix object x
  x1 <- as.data.frame(t(x))
  marker.vec <- row.names(x1)
  
  # First create a data.frame of metadata
  meta_df <- list(marker = marker.vec, chrom = chrom.vec, pos = pos.vec, cmPos = cmPos.vec, allele = allele.vec)
  meta_df <- as.data.frame(subset(meta_df, !sapply(meta_df, is.null)), stringsAsFactors = FALSE)
  
  x2 <- cbind(meta_df, x1)
  row.names(x2) <- NULL
  
  # Convert to marker_geno and return
  as_marker_geno(x = x2, chrom.col = "chrom", pos.col = "pos", cmPos.col = "cmPos", allele.col = "allele")
  
}





#'
#' @describeIn as_marker_geno
#' 
#' @param file Path to the CSV file. The first column must include the marker names; subsequent column names must be marker metadata (chromosom, position, etc.) or sample names. 
#' @param ... Additional arguments to pass to \code{\link[readr]{read_csv}}.
#' 
#' @importFrom readr read_csv
#' 
#' @export
#'
read_geno <- function(file, chrom.col, pos.col, cmPos.col, allele.col, ...) {
  
  # Error checking
  stopifnot(file.exists(file))
  stopifnot(is.character(chrom.col))
  stopifnot(is.character(pos.col))
  if (!missing(cmPos.col)) stopifnot(is.character(cmPos.col)) else cmPos.col <- NULL
  if (!missing(allele.col)) stopifnot(is.character(allele.col)) else allele.col <- NULL
  
  # Capture additional arguments
  other.args <- list(...)
  
  # Read in the file using read_csv
  data_in <- do.call(what = "read_csv", args = c(list(file = file), other.args))
  
  # Coerce to marker.geno object
  as_marker_geno(x = data_in, chrom.col = chrom.col, pos.col = pos.col, cmPos.col = cmPos.col,
                 allele.col = allele.col)
  
}


#' Manipulate marker.geno objects
#' 
#' 
#' @param ... Two or more \code{marker.geno} objects.
#' @param lst A list of \code{marker.geno} objects.
#' 
#' @import purrr
#' @importFrom dplyr full_join as_tibble
#' 
#' @export
#' 
#' 
merge_geno <- function(..., lst = list()) {
  
  ## Error checking
  stopifnot(is.list(lst))
  
  # Combine the dots with the list
  geno_list <- c(list(...), lst)
  
  # If only one object was passed, return it
  if (length(geno_list) == 1) return(geno_list[[1]])
  
  # Make sure each inherits the correct class
  if (!all(sapply(X = geno_list, FUN = inherits, "marker.geno"))) stop("One or more of the objects passed are not 'marker.geno' objects.")
  # Make sure none of the elements has duplicates
  if (!all(sapply(X = lapply(X = geno_list, FUN = attr, "duplicates"), is.null))) {
    stop("One or more of the objects passed has/have duplicate sample names; please resolve the duplicates before merging.")
  }
  
  # Merge based on common metadata columns
  col_list <- lapply(geno_list, attr, "cols")
  common_cols <- Reduce(f = intersect, x = col_list)
  
  # Full merge
  # merge1 <- function(x, y) merge.data.frame(x = x, y = y, by = common_cols, no.dups = TRUE, all = TRUE)
  merge1 <- function(x, y) full_join(x = x, y = y, by = common_cols)
  geno_merged <- Reduce(f = merge1, x = geno_list)
  
  ## Rename the duplicates
  # Remove the suffixes to get to the base sample names
  base_sample_names <- gsub(pattern = "(.y){1,}$|(.x){1,}$", replacement = "", x = names(geno_merged))
  
  # Make the names unique
  unique_base_sample_names <- make.unique(base_sample_names, sep = ";")
  # Add the names back to the data.frame
  names(geno_merged) <- unique_base_sample_names
  
  # Identify the duplicates; store this information as an attribute in the df
  dups <- base_sample_names[duplicated(base_sample_names)]
  dups <- if (length(dups) == 0) NULL else dups
  
  
  # Convert to tibble
  geno_merged1 <- as_tibble(geno_merged)
  
  # Add other attributes
  col_list1 <- col_list[[which.max(sapply(X = col_list, FUN = length))]]
  class(geno_merged1) <- union("marker.geno", class(geno_merged1))
  attr(geno_merged1, "cols") <- subset(col_list1, col_list1 %in% common_cols)
  attr(geno_merged1, "duplicates") <- dups
  
  return(geno_merged1)
  
}

#' Resolve duplicate entries in geno objects
#' 
#' @param x A \code{marker.geno} or \code{marker.matrix} object with duplicates.
#' 
#' @details 
#' Duplicates are resolved by using the consensus genotype (i.e. the most frequent genotype at a marker across duplicate samples). If more than 2 marker states are present among the duplicate samples, the genotype call is set to missing (i.e. NA).
#' 
#' @export
#' 
resolve_duplicate_geno <- function(x) {
  UseMethod("resolve_duplicate_geno", x)
}

#' 
#' @rdname resolve_duplicate_geno
#' 
#' @export
#' 
resolve_duplicate_geno.marker.geno <- function(x) {
  
  # If no duplicates, return the original data.frame
  dups <- attr(x, "duplicates")
  if (is.null(dups)) {
    warning("No duplicates were found; returning original marker.geno object.\n")
    return(x)
    
  }
  

  # Subset the duplicated data
  dup_cols <- unlist(lapply(X = lapply(X = dups, FUN = startsWith, x = names(x)), FUN = which))
  dup_data <- x[dup_cols]
  non_dup_data <- x[-dup_cols]
  
  # Create a list to store resolved duplicate data
  resolved_dups <- list()
  # Alternative rounding function
  round2 <- function(x) ifelse(x <= -0.5, -1, round(x))
  
  # Iterate over these names and merge the corresponding data
  for (dup_name in dups) {
    # Pull out corresponding data
    dup_name_data <- dup_data[startsWith(x = names(dup_data), prefix = dup_name)]
    
    # The consensus genotype is the most frequent, calculated as the mean
    mean_calls <- round2(rowMeans(x = dup_name_data, na.rm = TRUE))
    
    if (ncol(dup_name_data) > 2) {
      # Look for snps with > 2 unique calls; make these NA
      nDistinct <- apply(X = dup_name_data, MARGIN = 1, FUN = function(snp) length(unique(snp)))
      mean_calls[nDistinct > 2] <- NA
      
    }
    
    # Add the consensus call to the list
    resolved_dups[[dup_name]] <- mean_calls
    
  }
  
  # cbind
  resolved_dup_data <- as.data.frame(do.call("cbind", resolved_dups))
  # Merge with the no dup data
  geno_merged1 <- as_tibble(cbind(non_dup_data, resolved_dup_data))
  
  # Add other attributes
  class(geno_merged1) <- union("marker.geno", class(geno_merged1))
  attr(geno_merged1, "cols") <- attr(x, "cols")
  
  # return
  return(geno_merged1)
  
}

#' 
#' @rdname resolve_duplicate_geno
#' 
#' @export
#' 
resolve_duplicate_geno.marker.matrix <- function(x) {

  # Find duplicates
  sample_names <- rownames(x)
  # Make unique
  sample_names_unique <- make.unique(names = sample_names, sep = ";")
  
  # Identify the duplicates; store this information as an attribute in the df
  dups <- unique(gsub(pattern = ";[0-9]{1,}$", replacement = "", x = grep(pattern = ";[0-9]{1,}$", x = sample_names_unique, value = TRUE)))
  dups <- if (length(dups) == 0) NULL else dups
  
  if (is.null(dups)) {
    warning("No duplicates were found; returning original marker.geno object.\n")
    return(x)
    
  }
  
  
  # Subset the duplicated data
  dup_rows <- unlist(lapply(X = lapply(X = dups, FUN = startsWith, x = sample_names), FUN = which))
  dup_data <- x[dup_rows, ,drop = FALSE]
  non_dup_data <- x[-dup_rows, ,drop = FALSE]
  
  # Create a matrix to store the results
  resolved_dups <- matrix(data = NA, nrow = length(dups), ncol = ncol(x), dimnames = list(dups, colnames(x)))
  # Alternative rounding function
  round2 <- function(x) ifelse(x >= 0.5, 1, ifelse(x <= -0.5, -1, round(x)))
  roundHomo <- function(x) ifelse(x < 0, -1, 1)

  # Iterate over these names and merge the corresponding data
  for (dup_name in dups) {
    # Pull out corresponding data
    dup_name_data <- dup_data[startsWith(x = rownames(dup_data), prefix = dup_name), ,drop = FALSE]
    
    # The consensus genotype is the most frequent, calculated as the mean
    mean_calls <- colMeans(x = dup_name_data, na.rm = TRUE)
    # If there are no zeros, round to nearest homozygote
    zeros <- colSums(dup_name_data == 0, na.rm = TRUE) > 0
    mean_calls[zeros] <- round2(mean_calls[zeros])
    mean_calls[!zeros] <- roundHomo(mean_calls[!zeros])

    if (nrow(dup_name_data) > 3) {
      # Look for snps with > 2 unique calls; make these NA
      nDistinct <- apply(X = dup_name_data, MARGIN = 2, FUN = function(snp) length(unique(snp)))
      mean_calls[nDistinct > 3] <- NA
      
    }
    
    # Add the consensus call to the matrix
    resolved_dups[dup_name,] <- mean_calls
    
  }
  
  # cbind
  # Merge with the no dup data
  geno_merged1 <- rbind(non_dup_data, resolved_dups)
  
  # Add other attributes
  class(geno_merged1) <- union("marker.matrix", class(geno_merged1))
  # return
  return(geno_merged1)
  
}






#' @describeIn merge_geno
#' Edit marker genotype sample names
#'
#' @param replacement A \code{data.frame} with two columns: the sample name as specified in the \code{marker.geno} object; and the corrected # sample name to use.
#'
#' @export
#' 
edit_samples_geno <- function(x, replacement) {

  # Error checking
  stopifnot(inherits(x, "marker.geno"))
  stopifnot(is.data.frame(replacement))
  
  # Meta columns
  meta_col <- attr(x, "cols")
  
  # Subset the non-meta columns
  x1 <- x[! names(x) %in% meta_col]
  xMeta <- x[names(x) %in% meta_col]
  
  
  ## Make corrections
  # Find the index of names of x in the first replacement column
  ind <- match(x = names(x1), table = replacement[[1]])
  # Use this index to find the replacement names
  repl <- replacement[[2]][ind]
  # any NA replacement names become the original names
  repl[is.na(repl)] <- names(x1)[is.na(repl)]
  # Add these replacement names
  names(x1) <- repl
  
  # Recombine xMeta and the new x1
  x2 <- cbind(xMeta, x1)
  
  name.col <- meta_col[["name.col"]]
  chrom.col <- meta_col[["chrom.col"]]
  pos.col <- meta_col[["pos.col"]]
  cmPos.col <- meta_col["cmPos.col"]
  cmPos.col <- if (is.na(cmPos.col)) NULL else cmPos.col
  allele.col <- meta_col["allele.col"]
  allele.col <- if (is.na(allele.col)) NULL else allele.col
  
  as_marker_geno(x = x2, chrom.col = chrom.col, pos.col = pos.col, cmPos.col = cmPos.col, allele.col = allele.col)
  
}


#' Convert objects to marker.matrix objects
#' 
#' @param x A \code{marker.geno} object or a matrix of marker genotypes.
#' 
#' @return 
#' A \code{marker.matrix} object.
#' 
#' @export
#' 
as_marker_matrix <- function(x, ...) {
  UseMethod(generic = "as_marker_matrix", object = x)
}

#' 
#' @rdname as_marker_matrix
#' 
#' @export
#' 
as_marker_matrix.marker.geno <- function(x) {
  
  # Find the metadata columns
  meta_cols <- attr(x, "cols")
  
  # Get the name.col; remove others
  x1 <- x[! names(x) %in% meta_cols[! names(meta_cols) %in% "name.col"]]
  # Pull out the metadata information
  x1_meta <- x[names(x) %in% meta_cols]
  
  # Convert to matrix
  x2 <- as.data.frame(x1)
  row.names(x2) <- x2[[1]]
  x3 <- x2[-1]
  
  x4 <- t(as.matrix(x3))
  # Add class and attributes
  attr(x4, "snp.metadata") <- as.data.frame(x1_meta)
  class(x4) <- union("marker.matrix", class(x4))
  return(x4)
  
}

#' 
#' @rdname as_marker_matrix
#' 
#' @param chrom.vec A character vector of chromosome assignments; the length of this vector must equal the number of markers.
#' @param pos.vec A numeric vector of positions; the length of this vector must equal the number of markers.
#' @param cmPos.vec A numeric vector of genetic positions; the length of this vector must equal the number of markers.
#' @param allele.vec A character vector of marker alleles; the length of this vector must equal the number of markers.
#' 
#' @export
#' 
as_marker_matrix.matrix <- function(x, chrom.vec, pos.vec, cmPos.vec = NULL, allele.vec = NULL) {
  # Convert to marker geno
  x1 <- as_marker_geno.marker.matrix(x = x, chrom.vec = chrom.vec, pos.vec = pos.vec, cmPos.vec = cmPos.vec, allele.vec = allele.vec)
  # Convert to marker matrix and return
  as_marker_matrix.marker.geno(x = x1)
  
}



#' Summarize marker genotype data
#' 
#' @param x A \code{marker.matrix} object or \code{marker.geno} object.
#' @param simple Logical; should the summary be simple (i.e. just the number of markers and samples)? If FALSE, a full summary is calculated.
#' 
#' @details 
#' Calculates the number of markers, number of samples, marker and sample missingness,
#' marker and sample heterozygosity, marker allele frequency, and minor allele frequency.
#' 
#' @export
#' 
summary_geno <- function(x, simple = FALSE) {
  UseMethod("summary_geno", x)
}

#' 
#' @rdname summary_geno
#' 
#' @export
#' 
summary_geno.marker.matrix <- function(x, simple = FALSE) {

  stopifnot(is.logical(simple))
  
  # Basic summaries - number of markers and entries
  nMar <- ncol(x)
  nInd <- nrow(x)
  
  if (!simple) {
    # Missingness
    isNA <- is.na(x)
    marker_miss <- colMeans(isNA)
    ind_miss <- rowMeans(isNA)
    
    # MAF
    af <- colMeans(x + 1, na.rm = TRUE) / 2
    maf <- pmin(af, 1 - af)
    
    # Heterozyosity
    isHet <- x == 0
    marker_het <- colMeans(isHet, na.rm = TRUE)
    ind_het <- rowMeans(isHet, na.rm = TRUE)
    
    # ## Adjacent marker LD - calculated as r^2
    # # Pull out the snp metadata
    # snp_meta <- attr(x, "snp.metadata")
    # meta_cols <- attr(x, "cols")
    # # Order by chrom and pos; split by chrom
    # snp_meta1 <- snp_meta[order(snp_meta[[meta_cols[["chrom.col"]]]], snp_meta[[meta_cols[["pos.col"]]]]),]
    # snp_meta_split <- split(x = snp_meta1, f = snp_meta1$chrom)
    # # Get marker names
    # marker_names_split <- lapply(X = snp_meta_split, FUN = "[[", "marker")
    # 
    # # Create an empty list
    # marker_ld_list <- vector("list", length = length(marker_names_split))
    # # Iterate over chromosomes
    # for (i in seq_along(marker_names_split)) {
    #   markers_i <- marker_names_split[[i]]
    #   marker_mat_i <- x[,markers_i, drop = FALSE]
    #   
    #   # Create an empty matrix
    #   marker_ld_mat_i <- matrix(data = NA, nrow = ncol(marker_mat_i), ncol = ncol(marker_mat_i), 
    #                             dimnames = list(colnames(marker_mat_i), colnames(marker_mat_i)))
    #   
    #   # Iterate
    #   for (j in seq(2, ncol(marker_ld_mat_i))) {
    #     j1 <- j - 1
    #     # Calculate and store the correlation
    #     marker_ld_mat_i[j1, j] <- cor(marker_mat_i[,c(j1, j)], use = "complete.obs")[1,2]
    #     
    #   }
    #   
    #   # Add the matrix to the list
    #   # Calculate the squared correlation coefficient
    #   marker_ld_list[[i]] <- marker_ld_mat_i[,-1]^2
    #   
    # } # Close the per-chrom loop
    
  } else {
    marker_miss <- ind_miss <- marker_het <- ind_het <- af <- maf <- marker_ld_list <- NULL
    
  }
  
  
  # Assemble into a list
  summ_out <- list(
    nMarker = nMar,
    nSample = nInd,
    marker.missing = marker_miss,
    sample.missing = ind_miss,
    marker.het = marker_het,
    sample.het = ind_het,
    allele.freq = af,
    minor.allele.freq = maf
  )
  
  # Add two special classes
  class(summ_out) <- c("marker.summary", "marker.matrix.summary")
  return(summ_out)
  
}


#' 
#' @rdname summary_geno
#' 
#' @export
#' 
summary_geno.marker.geno <- function(x, simple = FALSE) {
  
  stopifnot(is.logical(simple))
  
  
  # Convert to marker matrix and calculate summary
  marker.mat <- as_marker_matrix.marker.geno(x)
  attr(marker.mat, "cols") <- attr(x, "cols")
  summ_out <- summary_geno.marker.matrix(x = marker.mat, simple = simple)
  
  # Replace the class type
  class(summ_out) <- c("marker.summary", "marker.geno.summary")
  return(summ_out)
  
}

#' 
#' @rdname summary_geno
#' 
#' @export
#' 
print.marker.summary <- function(x) {
  
  # Show
  cat("\nMarker summary \n\n")
  cat("Number of samples: ", x$nSample, "\n")
  cat("Number of markers: ", x$nMarker, "\n")
  
  # Test if simple
  if (!is.null(x$marker.missing)) {
  
    # Merge common summaries
    suppressWarnings({
      miss_summ <- rbind(
        markers = summary(x$marker.missing),
        samples = summary(x$sample.missing)
      )
    })
    
    suppressWarnings({
      het_summ <- rbind(
        markers = summary(x$marker.het),
        samples = summary(x$sample.het)
      )
    })
    
    cat("\nDistribution of marker/sample missingness:\n")
    print(miss_summ)
    
    cat("\nDistribution of marker/sample heterozygosity:\n")
    print(het_summ)
    
    cat("\nDistribution of marker minor allele frequency:\n")
    print(summary(x$minor.allele.freq))
    
  }
    
}


#' 
#' @rdname summary_geno
#' 
#' @export
#' 
print.marker.geno <- function(x) {
  cat("Printing a summary of this marker.geno object...\n\n")
  print(summary_geno(x, simple = TRUE))
}

#' 
#' @rdname summary_geno
#' 
#' @export
#' 
print.marker.matrix <- function(x) {
  cat("Printing a summary of this marker.matrix object...\n\n")
  print(summary_geno(x, simple = TRUE))
}




#' 
#' @rdname summary_geno
#' 
#' @export
#' 
plot.marker.summary <- function(x) {
  
  ## If the summary was simple (i.e. no x$marker.missing or other non-simple summaries), return invisible and warning
  if (is.null(x$marker.missing)) {
    # Print warning
    warning("The 'marker.summary' object does not contain plotable summaries.")
    
    # Return invisible
    invisible(NULL)
    
  }
  
  # Get the default par
  opar <- par()

  # Set par
  par(mfrow = c(3, 2))
  # Reset to not ask the user
  devAskNewPage(FALSE)

  # Plot
  hist(x$marker.missing, main = "Marker missingness", xlab = NULL)
  # Plot
  hist(x$sample.missing, main = "Sample missingness", xlab = NULL)
  # Plot
  hist(x$marker.het, main = "Marker heterozygosity", xlab = NULL)
  # Plot
  hist(x$sample.het, main = "Sample missingness", xlab = NULL)
  # Plot
  hist(x$minor.allele.freq, main = "Marker minor allele frequency", xlab = NULL)
  
  ## Adjacent marker LD

  
  adj_ld <- lapply(X = x$adjacent.ld, FUN = diag)
  
  # New plot - ask before resetting
  devAskNewPage(TRUE)
  plot.new()
  
  # Plot histograms
  for (i in seq_along(adj_ld)) {
    devAskNewPage(ask = i == 1)
    # Plot
    hist(x = adj_ld[[i]], main = paste0("Chromosome ", i, " adjacent marker LD"), xlab = NULL)
    
  }
  
  
  # Reset par
  suppressWarnings(par(opar))
  # Reset to false
  devAskNewPage(FALSE)
  
  # Return invisible
  invisible(NULL)
  
}


#' Filter marker genotype data
#' 
#' @param x A \code{marker.matrix} or \code{marker.geno} object.
#' @param marker.summary An optional marker summary object for x. If passed, it may speed up the filtering function.
#' @param max.marker.miss The maximum missingness rate for a marker.
#' @param max.sample.miss The maximum missingness rate for a sample.
#' @param range.marker.het A length-two vector of the minimum and maximum marker heterozygosity rate.
#' @param range.sample.het A length-two vector of the minimum and maximum sample heterozygosity rate.
#' @param min.maf The minimum minor allele frequency for a marker.
#' @param markers A character vector of markers to select. If passed, can be \code{"all"} or a character vector of markers to subset; in this case, all other filtering criteria are ignored.
#' @param samples A character vector of samples to select. If passed, can be \code{"all"} or a character vector of samples to subset; in this case, all other filtering criteria are ignored.
#' 
#' 
#' 
#' @export
#' 
filter_geno <- function(x, marker.summary, max.marker.miss = 1, max.sample.miss = 1, range.marker.het = c(0, 1),
                        range.sample.het = c(0, 1), min.maf = 0, markers, samples, verbose = TRUE) {
  
  stopifnot(inherits(x, c("marker.geno", "marker.matrix")))
  
  # Error checking
  stopifnot(is.numeric(max.marker.miss))
  stopifnot(is.numeric(max.sample.miss))
  stopifnot(is.numeric(range.marker.het))
  stopifnot(is.numeric(range.sample.het))
  stopifnot(is.numeric(min.maf))
  # stopifnot(is.numeric(max.adj.marker.ld))
  stopifnot(is.logical(verbose))
  
  if (!missing(markers)) {
    stopifnot(is.character(markers))
  } else {
    markers <- NULL
  }
  if (!missing(samples)) {
    stopifnot(is.character(samples))
  } else {
    samples <- NULL
  }
  
  if (!missing(marker.summary)) stopifnot(inherits(marker.summary, "marker.summary"))
  
  stopifnot(max.marker.miss >= 0 & max.marker.miss <= 1)
  stopifnot(max.sample.miss >= 0 & max.sample.miss <= 1)
  stopifnot(min.maf >= 0 & min.maf <= 1)
  # stopifnot(max.adj.marker.ld >= 0 & max.adj.marker.ld <= 1)
  
  stopifnot(length(range.marker.het) == 2)
  if (range.marker.het[1] > range.marker.het[2]) 
    stop("The first element of range.marker.het is the minimum heterozygosity; the second element is the maximum.")
  stopifnot(length(range.sample.het) == 2)
  if (range.sample.het[1] > range.sample.het[2]) 
    stop("The first element of range.sample.het is the minimum heterozygosity; the second element is the maximum.")
  
  # Do not summarize the object if markers and sample are BOTH passed
  if (is.null(markers) | is.null(samples)) {
    # First see if the marker summary was already passed
    if (!missing(marker.summary)) {
      x_summ <- marker.summary
    
    } else {
      x_summ <- summary_geno(x)
    }
    
  } 
  
  if (inherits(x, "marker.geno")) {
    x_markers <- x[[1]]
    x_samples <- setdiff(names(x), attr(x, "cols"))
    
  } else if (inherits(x, "marker.matrix")) {
    x_markers <- colnames(x)
    x_samples <- rownames(x)
    
  }
  
  # Count the number of markers and samples in each object 
  nMar <- length(x_markers)
  nSamp <- length(x_samples)
    
  # Apply the filters
  if (is.null(markers)) {
    markers_select_ind <- x_summ$marker.missing <= max.marker.miss & x_summ$marker.het >= range.marker.het[1] & 
      x_summ$marker.het <= range.marker.het[2] & x_summ$minor.allele.freq >= min.maf
    markers_select_ind <- which(markers_select_ind)
    
    # # Look at adjacent marker LD
    # adj_ld <- x_summ$adjacent.ld
    # # For each chromosome look for pairs of markers with elevated marker LD
    # marker_remove <- list()
    # for (i in seq_along(adj_ld)) {
    #   adj_ld_mat_i <- adj_ld[[i]]
    #   # Which pairs were above the threshold
    #   ld_above_max <- which(diag(adj_ld_mat_i) > max.adj.marker.ld)
    #   # Remove any that border
    #   ld_above_max1 <- subset(ld_above_max, !c(FALSE, diff(ld_above_max) == 1))
    #   
    #   # Get the second marker in each pair and add to the list
    #   marker_remove[[i]] <- colnames(adj_ld_mat_i)[ld_above_max1]
    #   
    # } # Close the loop
    # 
    # # Unlist marker_remove
    # marker_remove1 <- unlist(marker_remove)
    # # Subset markers_select_ind
    # markers_select_ind <- subset(markers_select_ind, !names(markers_select_ind) %in% marker_remove1)
    
    
  } else if (all(markers %in% "all")) {
    markers_select_ind <- seq_along(x_markers)
    
  } else {
    markers_select_ind <- which(x_markers %in% markers)
    
  }
  
  if (is.null(samples)) {
    sample_select_ind <- x_summ$sample.missing <= max.sample.miss & x_summ$sample.het >= range.sample.het[1] & 
      x_summ$sample.het <= range.sample.het[2]
    sample_select_ind <- which(sample_select_ind)
    
  } else if (all(samples %in% "all")) {
    sample_select_ind <- seq_along(x_samples)
    
  } else {
    sample_select_ind <- which(x_samples %in% samples)
    
  }
  
  ## Filter
  if (inherits(x, "marker.geno")) {
    x1 <- x[markers_select_ind, c(seq_along(attr(x, "cols")), sample_select_ind + length(attr(x, "cols"))), drop = FALSE]
    
  } else if (inherits(x, "marker.matrix")) {
    x1 <- x[sample_select_ind, markers_select_ind, drop = FALSE]
    
  }
  
  # Print notes about filtration, if called
  if (verbose) {
    cat(paste0("\nNumber of markers retained: ", length(markers_select_ind), " / ", nMar, ".\n"))
    cat(paste0("Number of samples retained: ", length(sample_select_ind), " / ", nSamp, ".\n"))
    
  }
  
  class(x1) <- class(x)
  
  # Return
  return(x1)
  
}


#' Prune SNPs based on linkage disequilibrium (LD)
#' 
#' @param x A \code{marker.geno} object.
#' @param maxLD The maximum LD of SNPs within a window.
#' @param posWindowSize The size of the sliding window (in bp).
#' @param snpWindowSize The maximum number of SNPs in the sliding window.
#' @param verbose Logical; should the function print status messages?
#' 
#' @return 
#' A \code{marker.geno} object with SNPs filtered according to LD.
#' 
#' @import SNPRelate
#' 
#' @export
#' 
prune_LD_geno <- function(x, maxLD = 0.95, posWindowSize, snpWindowSize, verbose = FALSE) {
  
  # Error checking
  stopifnot(inherits(x, "marker.geno"))
  stopifnot(is.numeric(maxLD))
  stopifnot(maxLD <= 1, maxLD >= 0)
  
  if (missing(posWindowSize) & missing(snpWindowSize)) {
    stop("You must specify a sliding window size by position ('posWindowSize') or number of SNPs ('snpWindowSize').")

  } else if (missing(posWindowSize)) {
    stopifnot(is.numeric(snpWindowSize), snpWindowSize > 0)
    posWindowSize <- NA
    
  } else if (missing(snpWindowSize)) {
    stopifnot(is.numeric(posWindowSize), posWindowSize > 0)
    snpWindowSize <- NA
    
  }
  
  stopifnot(is.logical(verbose))
  
  # Convert the marker geno object to a matrix
  x_mat <- as_marker_matrix.marker.geno(x = x)
  
  ## Convert the genotype data to a gds format
  # Vectors for metadata information
  snp_meta <- attr(x_mat, "snp.metadata")
  markers <- snp_meta[[1]]
  chroms <- as.numeric(as.factor(snp_meta[[2]]))
  pos <- snp_meta[[3]]
  
  #create gds file, can only do this 1 time per r session
  snpgdsCreateGeno(gds.fn = "full.gds", genmat = x_mat + 1, # Convert to 0, 1, 2 format
                   snp.id = markers, snp.chromosome = chroms, snp.position = pos,
                   snpfirstdim = FALSE)
  
  # Read the GDS back in
  full_gds_obj <- snpgdsOpen("full.gds", readonly = FALSE, allow.duplicate = FALSE, allow.fork = FALSE)
        
  # Implement the pruning - start at first
  if (verbose) cat("\nPruning LD from the beginning of chromosomes...\n\n")
  snpset_begin <- snpgdsLDpruning(gdsobj = full_gds_obj, remove.monosnp = TRUE, method = "corr", 
                                  slide.max.bp = posWindowSize, slide.max.n = snpWindowSize, ld.threshold = maxLD, 
                                  verbose = verbose, start.pos = "first")
  
  # Next go from the last
  if (verbose) cat("\nPruning LD from the ending of chromosomes...\n\n")
  snpset_end <- snpgdsLDpruning(gdsobj = full_gds_obj, remove.monosnp = TRUE, method = "corr", 
                                slide.max.bp = posWindowSize, slide.max.n = snpWindowSize, ld.threshold = maxLD, 
                                verbose = verbose, start.pos = "last")
  
  # Find the intersection of the SNP sets
  snpset_use <- intersect(unlist(snpset_begin), unlist(snpset_end))
  
  # Close the gds file and remove it
  snpgdsClose(gdsobj = full_gds_obj)
  invisible(file.remove("full.gds"))
  
  # Subset the marker geno file and return in the class in which it was submitted
  if (verbose) cat("\nFiltering the 'marker.geno' object and returning...\n")
  x_mat2 <- filter_geno(x = x, markers = snpset_use, samples = row.names(x_mat), verbose = verbose)
  
  # Return
  return(x_mat2)
  
}
