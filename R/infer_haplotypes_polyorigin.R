#' Reconstruct haplotypes using PolyOrigin
#'
#' @param geno A data.frame of genotype calls, where the first three columns are "marker," "chrom,", and "cM" (genetic position);
#' optionally, the fourth column can be the physical position with column name "bp". Subsequent columns are genotype calls for individuals.
#' @param ped A pedigree data frame where the first four columns are the individual name, the family (or population) name,
#' the mother genotype ID, and the father genotype ID.
#' @param out.prefix The prefix (i.e. no file extension) for saving the output of PolyOrigin.
#' @param ploidy Ploidy level; must be 2, 4, or 6.
#' @param filter.mars Logical; should the data for each family be filtered using mappoly?
#' @param max.ind.miss The maximum rate of missingness to retain an individual
#' @param max.mar.miss The maximum rate of missingness to retain an marker
#' @param julia.path The path to the Julia executable on the system.
#'
#' @importFrom mappoly read_geno_csv filter_missing
#' @import PolyOriginR
#'
#' @export
#'
infer_haplotypes_polyorigin <- function(geno, ped, out.prefix, ploidy = 2, filter.mars = FALSE, max.ind.miss = 0.20, max.mar.miss = 0.20,
                                        julia.path = "") {

  require(mappoly)

  # Error handling
  stopifnot(is.data.frame(geno))
  cols.match <- names(geno)[1:3] == c("marker", "chrom", "position")
  if (any(!cols.match)) {
    stop("The first 3 columns of 'geno' should be 'marker', 'chrom', and 'position'.")
  }

  stopifnot(is.data.frame(geno))
  stopifnot(is.data.frame(ped))
  stopifnot(is.character(out.prefix))
  stopifnot(is.logical(filter.mars))

  # Write a temporary folder
  if (!dir.exists("tmp")) dir.create("tmp")

  stopifnot(all(ploidy %in% c(2, 4, 6)))
  if (length(ploidy) > 1) {
    stopifnot(length(ploidy) != nrow(ped))
  }

  # Prepare the pedigree information
  ped1 <- ped[,c("id", "pop", "mother", "father")]
  ped1$ploidy <- ploidy

  ped1$mother <- ifelse(is.na(ped1$mother), "0", ped1$mother)
  ped1$father <- ifelse(is.na(ped1$father), "0", ped1$father)

  # If the geno information contains bp, save it for later
  if (names(geno)[4] == "bp") {
    map1 <- geno[,1:4]
    geno1 <- geno[,-4]
  } else {
    map1 <- geno[,1:3]
    geno1 <- geno
  }

  # Remove any markers with NA positions
  idx_keep <- which(!is.na(geno1[,3]))
  geno2 <- geno1[idx_keep, ]
  map2 <- map1[idx_keep, ]

  # Rename
  names(geno2)[1:3] <- c("marker", "chromosome", "pos")

  # Intersect genotypes in the pedigree with those in the genotype matrix
  common_genos <- intersect(ped1$id, names(geno2))

  geno3 <- geno2[,c("marker", "chromosome", "pos", common_genos)]
  ped2 <- ped1[ped1$id %in% common_genos, ]
  # ped2 <- ped1
  names(ped2) <- c("individual", "population", "motherid", "fatherid", "ploidy")

  if (filter.mars) {

    # Iterate over the families and load the data into mappoly for filtering
    ped3 <- ped2[ped2$population != 0, ]
    ped3_split <- split(ped3, ped3$population)
    geno_filt_list <- list()
    for (i in seq_along(ped3_split)) {
      ped_i <- ped3_split[[i]]
      offspring <- ped_i$individual
      parents <- unique(unlist(ped_i[,c("motherid", "fatherid")]))
      # Format the genotypes for mappoly2
      geno3_i <- geno3
      geno3_i$ref <- geno3_i$alt <- as.character(NA)
      geno3_i <- geno3_i[,c("marker", parents, "chromosome", "pos", "ref", "alt", offspring)]
      # Rename
      geno_names <- c(parents, offspring)
      geno_names <- gsub(pattern = "-", replacement = "_", x = geno_names)
      names(geno3_i) <- gsub(pattern = "-", replacement = "_", x = names(geno3_i))
      filename_i <- file.path("tmp", paste0(paste0(parents, collapse = "x"), ".csv"))
      write_csv(x = geno3_i, file = filename_i, progress = FALSE)
      # Read in to mappoly
      mappoly_in <- read_geno_csv(file.in = filename_i, ploidy = 2, filter.non.conforming = TRUE, elim.redundant = TRUE, verbose = FALSE)
      mappoly_in <- filter_missing(input.data = mappoly_in, type = "marker", filter.thres = max.mar.miss, inter = FALSE)
      mappoly_in <- filter_missing(input.data = mappoly_in, type = "individual", filter.thres = max.ind.miss, inter = FALSE)
      # Determine the markers to remove
      markers_keep <- mappoly_in$kept
      indiv_keep <- c(geno_names[1:2], mappoly_in$ind.names)
      geno_filt_list[[i]] <- list(mar = markers_keep, ind = indiv_keep)
    }

    markers_keep <- Reduce(union, lapply(X = geno_filt_list, FUN = "[[", "mar"))
    ind_keep <- Reduce(union, lapply(X = geno_filt_list, FUN = "[[", "ind"))

    geno4 <- geno3
    geno4 <- geno4[geno4$marker %in% markers_keep, ]
    idx <- c(1:3, sort(match(ind_keep, gsub("-", "_", names(geno4)))))
    geno4 <- geno4[,idx]

  } else {
    markers_keep <- geno3$marker
    ind_keep <- ped2$individual

    geno4 <- geno3
    geno4 <- geno4[geno4$marker %in% markers_keep, ]
    idx <- c(1:3, sort(match(ind_keep, names(geno4))))
    geno4 <- geno4[,idx]

  }

  # Adjust the ped and map accordingly
  ped4 <- ped2[ped2$individual %in% names(geno4), ]

  # Adjust population names
  pop_fct <- apply(ped4[,c("motherid", "fatherid")], MARGIN = 1, FUN = paste0, collapse = "x")
  pop_fct <- factor(pop_fct, levels = unique(pop_fct))
  pop_fct_num <- as.numeric(pop_fct) - 1
  ped4$population <- ifelse(pop_fct_num == 0, "0", paste0("pop", pop_fct_num))


  ## Write temporary files
  write.csv(x = ped4, file = file.path("tmp/pedfile.csv"), quote = FALSE, row.names = FALSE)
  write.csv(x = geno4, file = file.path("tmp/genofile.csv"), quote = FALSE, row.names = FALSE)

  # # Run polyorigin
  # ans <- polyOriginR(genofile = file.path("tmp/genofile.csv"), pedfile = file.path("tmp/pedfile.csv"), juliapath = julia.path)

  # Create a julia script to run polyOrigin
  con <- file("tmp/po.jl", open = "write")
  writeLines("using PolyOrigin;", con)
  writeLines(sub("X", "tmp/genofile.csv", "genofile=\"X\";"), con)
  writeLines(sub("X", "tmp/pedfile.csv", "pedfile=\"X\";"), con)
  writeLines("polyOrigin(genofile,pedfile,isphysmap=false,refineorder=false,refinemap=false,outstem=\"tmp/imputed\");", con)
  close(con)
  system(paste0(dirname(julia.path), "/julia -t auto tmp/po.jl"))

  # Read in the output CSV
  geno_prob_df <- read.csv(file = "tmp/imputed_genoprob.csv", header = TRUE, as.is = TRUE, check.names = FALSE)
  # Read in the pedfile
  ped_df <- read.csv(file = "tmp/pedfile.csv")

  # Rename columns in both dfs
  marker_order <- geno_prob_df$marker
  idx <- match(x = marker_order, table = map1$marker)
  geno_prob_df1 <- cbind(map1[idx,], geno_prob_df[,-1:-3])

  ped_df1 <- ped_df[ped_df$population != 0,]
  ped_df1 <- ped_df1[,c("individual", "motherid", "fatherid")]
  names(ped_df1) <- c("id", "parent1", "parent2")

  # Save the files
  write.csv(x = geno_prob_df1, file = paste0(out.prefix, "_genofile.csv"), quote = FALSE, row.names = FALSE)
  write.csv(x = ped_df1, file = paste0(out.prefix, "_pedfile.csv"), quote = FALSE, row.names = FALSE)

  cat("\n\nHaplotype reconstruction complete.")
  cat("\nOutput genotype probabilities and pedigree file are saved at", out.prefix)

  # Delete the temporary directory
  unlink("tmp/", force = TRUE, recursive = TRUE)

}



#' Convert the output of polyorigin to a genoprobs object for R/qtl2
#'
#' @param genofile Path to the _genofile.csv output from polyorigin haplotype reconstruction
#' @param pedfile Path to the _pedfile.csv output from polyorigin haplotype reconstruction
#' @param ploidy The ploidy level
#'
#'
polyorigin_to_rqtl2 <- function(genofile, pedfile, ploidy = 2) {

  # error
  stopifnot(file.exists(genofile))
  stopifnot(file.exists(pedfile))

  # Read in the files
  geno_probs <- read.csv(file = genofile, header = TRUE, as.is = TRUE, check.names = FALSE)
  ped <- read.csv(file = pedfile, header = TRUE, as.is = TRUE, check.names = FALSE)

  Pars <- unique(ped[,c("parent1", "parent2")])
  par.hapl <- apply(X = Pars, MARGIN = 1, FUN = function(pars) {
    fam_pars <- paste0(rep(pars, each = 2), ".", c(1,2))
    c("1" = paste0(fam_pars[c(1,3)], collapse = ":"), "2" = paste0(fam_pars[c(2,3)], collapse = ":"),
      "3" = paste0(fam_pars[c(1,4)], collapse = ":"), "4" = paste0(fam_pars[c(2,4)], collapse = ":"))
  }, simplify = FALSE)

  par.hapl <- unique(unlist(par.hapl))

  # Create an array of genotype probabilities
  Ind <- ped$id
  nInd <- length(Ind)
  nMar <- nrow(geno_probs)
  MarChrom <- split(geno_probs$marker, geno_probs$chrom)
  nChrom <- length(MarChrom)
  nMarChrom <- sapply(MarChrom, length)

  geno_probs_array <- lapply(X = MarChrom, FUN = function(x) array(data = as.numeric(0), dim = c(nInd, length(par.hapl), length(x)), dimnames = list(Ind, par.hapl, x)))

  # Split the pedigree by family
  ped$pop <- as.numeric(as.factor(paste0(ped$parent1, ".", ped$parent2)))
  ped_split <- split(ped, ped$pop)

  # Iterate over the split pedigree
  for (fam in ped_split) {
    # The F1 codes are 1,2,3,4 corresponding to 1-3, 1-4, 2-3, 2-4 haplotype combinations
    #
    # First determine the haplotype combination
    fam_pars <- unique(unlist(fam[,c("parent1", "parent2")]))
    fam_pars <- paste0(rep(fam_pars, each = 2), ".", c(1,2))
    fam_haplo_code <- c("1" = paste0(fam_pars[c(1,3)], collapse = ":"), "2" = paste0(fam_pars[c(2,3)], collapse = ":"),
                        "3" = paste0(fam_pars[c(1,4)], collapse = ":"), "4" = paste0(fam_pars[c(2,4)], collapse = ":"))

    # Get the geno probs for these individuals
    geno_probs_fam <- geno_probs[,c("marker", "chrom", fam$id)]
    # Split by chromosome
    geno_probs_fam_split <- split(geno_probs_fam, geno_probs_fam$chrom)

    # Iterate over chromosomes
    for (i in seq_along(geno_probs_fam_split)) {
      chrom_i <- names(geno_probs_fam_split)[i]
      probs_i <- geno_probs_fam_split[[i]]
      mars_i <- probs_i$marker
      # Convert to a matrix
      probs_i <- as.matrix(probs_i[,-1:-2])
      row.names(probs_i) <- mars_i

      # Iterate over columns (individuals)
      for (j in seq_len(ncol(probs_i))) {
        ind <- probs_i[,j]
        indname <- colnames(probs_i)[j]
        ind_parse <- strsplit(x = ind, split = "=>", fixed = TRUE)
        for (k in seq_along(ind_parse)) {
          mar_k <- names(ind_parse)[k]
          x_parse <- strsplit(x = ind_parse[[k]], split = "|", fixed = TRUE)
          par_haplo_match <- fam_haplo_code[as.numeric(x_parse[[1]])]
          par_haplo_match <- match(x = par_haplo_match, table = par.hapl)

          # Add probabilities to the array
          geno_probs_array[[chrom_i]][indname, par_haplo_match, mar_k] <- as.numeric(x_parse[[2]])

        }
      }
    }
  }

  # Add attributes
  class(geno_probs_array) <- c("calc_genoprob", "list")
  attr(geno_probs_array, "crosstype") <- ""
  attr(geno_probs_array, "is_x_chr") <- setNames(rep(FALSE, length(geno_probs_array)), names(geno_probs_array))
  attr(geno_probs_array, "alleles") <- unique(unlist(strsplit(x = dimnames(geno_probs_array[[1]])[[2]], split = ":", fixed = TRUE)))
  attr(geno_probs_array, "alleleprobs") <- FALSE

  # Return the probabilities
  return(geno_probs_array)

}





