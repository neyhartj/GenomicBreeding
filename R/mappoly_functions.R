
# This function doesn't work yet
infer_haplotypes_mappoly <- function(geno, ped, ploidy = c(2, 4, 6), remove.uninformative = TRUE) {

  require(mappoly)

  stopifnot(is.data.frame(geno))
  stopifnot(is.data.frame(ped))
  ploidy <- match.arg(ploidy)

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

  # Prepare the pedigree information
  ped1 <- ped[,c("id", "pop", "mother", "father")]
  ped1$ploidy <- ploidy

  ped1$mother <- ifelse(is.na(ped1$mother), "0", ped1$mother)
  ped1$father <- ifelse(is.na(ped1$father), "0", ped1$father)

  # Split the pedigree by population
  ped2 <- ped1[ped1$pop != 0,]
  ped2_split <- split(ped2, ped2$pop)

  # Write temporary files
  if (!dir.exists("tmp")) dir.create("tmp")
  # Iterate over the populations
  for (ped_i in ped2_split) {
    # Character vector of parents
    parents <- unique(unlist(ped_i[,c("mother", "father")]))
    # Same for offspring
    offspring <- ped_i$id

    # Create a geno object
    geno_i <- geno2[,c("marker", parents, "chrom", "cM", offspring)]

    # Save the object
    tmp_geno_file <- file.path("tmp/genofile.csv")
    write.csv(x = geno_i, file = tmp_geno_file, quote = FALSE, row.names = FALSE)

    # Read the file into mappoly
    dat_i <- read_geno_csv(file.in = tmp_geno_file, ploidy = ploidy, elim.redundant = FALSE, verbose = FALSE)

    # Make a sequence object with the order of the markers preserved
    if (remove.uninformative) {
      keep <- dat_i$mrk.names
    } else {
      keep <- geno_i$marker
    }

    # Subset the map
    map_i <- geno_i[geno_i$marker %in% keep, c(1, 4, 5)]
    # Split the map into chromosomes
    map_i_split <- split(map_i, map_i$chrom)

    # Create a list of seq objects based on the order in the map
    seq_list <- lapply(map_i_split, function(df) {
      make_seq_mappoly(dat_i, df$marker)              # keeps order verbatim
    })

    # Calculate recombination fractions between markers using map data
    rf_list <- lapply(map_i_split, function(df) {
      if (nrow(df) < 2) return(numeric(0))
      mappoly:::mf_k(diff(df$cM))
    })

    # Create artificial mappoly map objects
    mappoly_maps_list <- mapply(function(seq.obj, df, rf.vec) {
      idx <- seq.obj$seq.num              # indices in the original data

      map.info <- list(
        ploidy      = dat_i$ploidy,
        n.mrk       = length(idx),
        seq.num     = idx,
        mrk.names   = df$marker,
        seq.dose.p1 = dat_i$dosage.p1[idx],
        seq.dose.p2 = dat_i$dosage.p2[idx],
        chrom       = df$chrom[1],
        genome.pos  = df$cM,
        data.name   = deparse(substitute(dat_i)),
        ph.thres    = 0
      )

      core <- list(seq.num = idx,
                   seq.rf  = rf.vec,
                   seq.ph  = NULL,
                   loglike = NA)

      structure(list(info = map.info, maps = list(core)),
                class = "mappoly.map")

    }, seq_list, map_i_split, rf_list, SIMPLIFY = FALSE)

    # Estimate linkage phases
    mappoly_maps_list1 <- lapply(mappoly_maps_list, est_full_hmm_with_global_error, error = 0.01, tol = 1e-3)

    test <- est_full_hmm_with_global_error(input.map = mappoly_maps_list$Vmac_chr01, error = 0.02, verbose = TRUE)


    # This preserves all markers, including those that may have been removed as "uninformative"
    seq_obj <- make_seq_mappoly(input.obj = dat_i, arg = keep)


    # Filter the map
    keep <- intersect(map_df$marker, dat$mrk.names)
    dat  <- drop_marker(dat, setdiff(dat$mrk.names, keep))
    map_df <- map_df[match(keep, map_df$marker), ]

  }


}





#' Convert haplotype probabilities from mappoly2
#'
#' @description
#' Functions to convert haplotype probabilities from mappoly2 for input into other
#' packages. Support is included for \code{diaQTL} and \code{qtl2}. The code to convert
#' to \code{diaQTL} format is borrowed from the \code{\link[diaQTL]{convert_mappoly}} function.
#'
#' @rdname export_mappoly2
#'
#' @param data An object of class \code{mappoly2.sequence} (for a single population) or \code{mappoly2.consensus.map}
#' (for multiple populations). The objects should be the output of \code{calc_haplotypes} or \code{calc_consensus_haplo}.
#' @param ploidy The ploidy level
#' @param out.prefix The prefix of the output files, including the path to the directory where the files should be saved
#'
#' @return
#' \code{export_mappoly2_to_diaqtl} - Creates a \code{genofile} and \code{pedfile} for import into the \code{diaQTL} package using \code{\link[diaQTL]{read_data}}.
#'
#' @examples
#'
#' \dontrun{
#'
#' alfalfa.f1 <- alfa_f1
#' alfalfa.f1 <- filter_data(alfalfa.f1, mrk.thresh = 0.2, ind.thresh = 0.1, plot.screening = FALSE)
#' alfalfa.f1 <- filter_individuals(alfalfa.f1, ind.to.remove = "F1.85.133", inter = FALSE, verbose = FALSE)
#' alfalfa.f1.all <- pairwise_rf(alfalfa.f1, mrk.scope = "all", ncpus = 8)
#' alfalfa.f1.all <- rf_filter(alfalfa.f1.all, thresh.LOD.ph = 5, thresh.LOD.rf = 5, thresh.rf = 0.15, probs = c(0.025, 0.975), diagnostic.plot = FALSE)
#' g <- group(x = alfalfa.f1.all, expected.groups = 8, comp.mat = TRUE, inter = FALSE)
#' s <- make_sequence(g, ch = list(1, 2, 3, 4, 5, 6, 7, 8))
#' s <- order_sequence(s, type = "genome")
#' s <- rf_filter(s, type = "genome", probs = c(0.025, 0.975), diag.markers = 50)
#' s <- pairwise_phasing(s, type = "genome", thresh.LOD.ph = 3, thresh.LOD.rf = 3, thresh.rf = 0.5, max.search.expansion.p1 = 10, max.search.expansion.p2 = 10)
#' s <- mapping(s, type = "genome", parent = "p1", ncpus = 8)
#' s <- mapping(s, type = "genome", parent = "p2", ncpus = 8)
#' s <- merge_single_parent_maps(s, type = "genome", ncpus = 8, error = 0.05)
#'
#' s <- calc_haplotypes(s, type = "genome", ncpus = 8)
#' s <- augment_phased_map(s, type = "genome", ncpus = 8)
#' s <- calc_haplotypes(s, type = "genome", ncpus = 8)
#' I195xJ432_map <- s
#'
#' # Load the second map
#' download.file("https://github.com/mmollina/mappoly2_vignettes/raw/main/I195_x_F1-85-209_map.rda",
#' destfile = "temp_file.rda")
#'
#' load("temp_file.rda")
#'
#' I195xF1_85_209_map <- mapping(I195xF1_85_209_map, type = "genome", error = 0.05, ncpus = 8)
#'
#' # Integrate the two maps
#' MAPs <- list(I195xJ432 = I195xJ432_map, I195xF1_85_209 = I195xF1_85_209_map)
#' prep.maps <- prepare_to_integrate(MAPs)
#'
#' consensus_map <- estimate_consensus_map(prep.maps, ncpus = 8, err = 0.05)
#' consensus_map <- calc_consensus_haplo(consensus_map, ncpus = 8)
#'
#'
#'
#' # Export to diaQTL
#' export_mappoly2_to_diaqtl(data = I195xJ432_map, ploidy = 4, out.prefix = "alfalfa_gp_")
#' export_mappoly2_to_diaqtl(data = consensus_map, ploidy = 4, out.prefix = "alfalfa_gp_consensus_")
#'
#' # Read into diaQTL
#' dat1 <- diaQTL::read_data(genofile = "alfalfa_gp_diaQTL_geno.csv", ploidy = 4, pedfile = "alfalfa_gp_diaQTL_ped.csv")
#' dat2 <- diaQTL::read_data(genofile = "alfalfa_gp_consensus_diaQTL_geno.csv", ploidy = 4, pedfile = "alfalfa_gp_consensus_diaQTL_ped.csv")
#'
#' # Export to qtl2
#' qtl2_dat1 <- export_mappoly2_to_qtl2(data = I195xJ432_map, ploidy = 4, type = "genome")
#' qtl2_dat2 <- export_mappoly2_to_qtl2(data = consensus_map, ploidy = 4, type = "genome")
#'
#' # Use in qtl2
#' kin1 <- qtl2::calc_kinship(probs = qtl2_dat2$probs, type = "loco")
#' pheno <- c(MAPs$I195xJ432$data$geno.dose['Alf_2262657_SNP1',] * 2,
#'            MAPs$I195xF1_85_209$data$geno.dose['Alf_2262657_SNP1',] * 2)
#' pheno <- as.matrix(pheno)
#'
#' ans <- qtl2::scan1(genoprobs = qtl2_dat2$probs, pheno = pheno, kinship = kin1, addcovar = qtl2_dat2$fam.covar)
#' plot(ans, map = qtl2_dat2$pmap)
#'
#' }
#'
#' @export
#'
export_mappoly2_to_diaqtl <- function(data, ploidy, type = c("genome", "mds"), out.prefix = "") {

  # Check class
  stopifnot(inherits(data, c("mappoly2.consensus.map", "mappoly2.sequence")))
  stopifnot(ploidy %in% c(2, 4))
  type <- match.arg(type)

  # Get the states
  if (ploidy == 4) {
    states <- c('1256', '1257', '1258', '1267', '1268', '1278', '1356', '1357',
                '1358', '1367', '1368', '1378', '1456', '1457', '1458', '1467',
                '1468', '1478', '2356', '2357', '2358', '2367', '2368', '2378',
                '2456', '2457', '2458', '2467', '2468', '2478', '3456', '3457',
                '3458', '3467', '3468', '3478')

  } else {
    states <- c("13", "14", "23", "24")
  }

  # Get the probabilities - a list of length nChr
  if (class(data) == "mappoly2.consensus.map") {
    prob_list <- lapply(data$consensus.map, FUN = "[[", "haploprob")
    # Create the pedigree data.frame
    ped <- data$consensus.map[[1]]$ph$pedigree
    # Parent names and numbers
    parents <- names(data$consensus.map[[1]]$ph$PH)
    parents <- setNames(seq_along(parents), parents)

    ped <- data.frame(id = row.names(ped), parent1 = names(parents[ped$Par1]), parent2 = names(parents[ped$Par2]),
                      row.names = NULL)

  } else {
    prob_list <- lapply(X = data$maps, FUN = function(x) x[[type]]$p1p2$hmm.phase[[1]]$haploprob)

    ped <- data.frame(id = data$data$ind.names, parent1 = data$data$name.p1, parent2 = data$data$name.p2, row.names = NULL)

  }

  outputAll <- NULL
  for (i in seq_along(prob_list)) {
    tmp <- as.matrix(prob_list[[i]])
    # Get the marker names for this lg
    if (class(data) == "mappoly2.consensus.map") {
      mars <- row.names(data$consensus.map[[i]]$ph$G)
      colnames(tmp)[-1:-3] <- mars
      ids <- as.character(row.names(tmp))
      ids <- unique(sub(pattern = "\\.Par[0-9]", replacement = "", x = ids))
    } else {
      mars <- row.names(data$maps[[i]][[type]]$p1p2$hmm.phase[[1]]$p1)
      colnames(tmp) <- c("parent", "ind", "homolog", mars)
      # Get id names
      ids <- data$data$screened.data$ind.names
      # Create rownames for tmp
      tmp_rownames <- paste0(ids[tmp[,"ind"]], ".Par", tmp[,"parent"])
      row.names(tmp) <- tmp_rownames
    }

    # Compute genotype probabilities from the haplotype probabilities
    states_idx <- lapply(strsplit(states, ""), as.numeric)
    names(states_idx) <- seq_along(states_idx)
    states_num <- sort(unique(unlist(states_idx)))
    state_indices_list <- lapply(states_num, function(x) seq(x, nrow(tmp), by = max(states_num)))

    tmp1 <- tmp[,-c(1:3)]
    geno_prob_states <- lapply(X = states_idx, FUN = function(ii) {
      states_ii <- state_indices_list[ii]
      Reduce(`*`, lapply(X = states_ii, function(iii) tmp1[iii, , drop = FALSE]))
    })
    geno_prob_states <- lapply(geno_prob_states, round, 3)
    geno_prob_states <- simplify2array(geno_prob_states)

    # Convert to data frame
    geno_prob_df <- as.data.frame.table(geno_prob_states)
    # Spread out the states
    geno_prob_df <- reshape(data = geno_prob_df, idvar = c("Var1", "Var2"), timevar = "Var3", direction = "wide")
    # Rename
    colnames(geno_prob_df) <- c("id", "marker", names(states_idx))
    # Remove zero probabilities
    prob <- apply(geno_prob_df[, -(1:2)], 1, function(x) x[x > 0] )
    geno_prob_df$prob <- sapply(prob, function(x) paste(paste0(names(x), collapse = "|"), paste0(x, collapse = "|"), sep = "=>") )

    # Spread ids out
    geno_prob_df <- geno_prob_df[,c("marker", "id", "prob")]
    geno_prob_df$id <- ids
    geno_prob_df <- reshape(data = geno_prob_df, idvar = "marker", timevar = "id", direction = "wide")
    names(geno_prob_df) <- gsub(pattern = "prob\\.", replacement = "", x = names(geno_prob_df))

    # Add chrom, bp, and cM names
    if (class(data) == "mappoly2.consensus.map") {
      chrom <- names(data$consensus.map)[i]
      cM <- cumsum(c(0, imf_h(data$consensus.map[[i]]$rf)))
      bp_list <- lapply(X = data$individual.maps, FUN = function(x) {
        bps <- x$maps[[i]]$genome$order
        data.frame(marker = row.names(bps), bps, row.names = NULL)
      })
      bp <- do.call(rbind, bp_list)
      row.names(bp) <- NULL
      bp <- unique(bp)

    } else {
      chrom <- unique(data$data$chrom[mars])
      cM <- cumsum(c(0, imf_h(data$maps[[i]][[type]]$p1p2$hmm.phase[[1]]$rf)))
      bp <- data$data$genome.pos[mars]
      bp <- data.frame(marker = mars, chrom = chrom, bp = bp, row.names = NULL)

    }

    output <- bp[match(x = unique(geno_prob_df$marker), bp$marker), ]
    output$cM <- cM
    names(output) <- c("marker", "chrom", "bp", "cM")
    output <- output[, c("marker", "chrom", "cM", "bp")]
    output <- data.frame(output, geno_prob_df[,-1])

    outputAll = rbind(outputAll, output)
  }

  write.csv(outputAll, file = paste0(out.prefix, "diaQTL_geno.csv"),  row.names = F)
  write.csv(ped, file = paste0(out.prefix, "diaQTL_ped.csv"), row.names = F)

}






#'
#' @rdname export_mappoly2
#'
#' @return
#' \code{export_mappoly2_to_qtl2} - A list with four elements:
#'    \itemize{
#'      \item{\code{probs} - An object of class \code{calc_genoprob} with allele probabilities for
#'      use in the \code{\link[qtl2]{qtl2}} package.}
#'      \item{\code{gmap} - A genetic map in list form for use in \code{\link[qtl2]{qtl2}}.}
#'      \item{\code{pmap} - A physical map in list form for use in \code{\link[qtl2]{qtl2}}.}
#'      \item{\code{fam.covar} - An incidence matrix for coding family effects when
#'      conducting a QTL scan in \code{\link[qtl2]{qtl2}}.}
#'  }
#'
#'
#' @export
#'
export_mappoly2_to_qtl2 <- function(data, ploidy, type = c("genome", "mds")) {

  # Check class
  stopifnot(inherits(data, c("mappoly2.consensus.map", "mappoly2.sequence")))
  stopifnot(ploidy %in% c(2, 4))
  type <- match.arg(type)

  # Get the probabilities - a list of length nChr
  if (class(data) == "mappoly2.consensus.map") {
    prob_list <- lapply(data$consensus.map, FUN = "[[", "haploprob")
    # Create the pedigree data.frame
    ped <- data$consensus.map[[1]]$ph$pedigree
    # Parent names and numbers
    parents <- names(data$consensus.map[[1]]$ph$PH)
    parents <- setNames(seq_along(parents), parents)

    ped <- data.frame(id = row.names(ped), parent1 = names(parents[ped$Par1]), parent2 = names(parents[ped$Par2]),
                      row.names = NULL)

    # Get the genetic and physical maps
    gmap <- lapply(X = data$consensus.map, FUN = function(x) {
      cm <- c(0, cumsum(imf_h(x$rf)))
      setNames(cm, row.names(x$ph$G))
    })
    names(gmap) <- seq_along(gmap)
    attr(gmap, "is_x_chr") <- setNames(rep(FALSE, length(gmap)), names(gmap))

    # Get the physical map
    genome_pos <- lapply(X = data$individual.maps, FUN = function(x) {
      pos <- x$data$genome.pos
      data.frame(marker = names(pos), bp = pos, row.names = NULL)
    })
    genome_pos <- unique(`row.names<-`(do.call(rbind, genome_pos), NULL))

    # Resolve issues where the same marker has different genome positions; take the first position
    genome_pos <- split(genome_pos, genome_pos$marker)
    genome_pos <- lapply(X = genome_pos, FUN = function(x) x[1,, drop = FALSE])
    genome_pos <- `row.names<-`(do.call(rbind, genome_pos), NULL)

    row.names(genome_pos) <- genome_pos[[1]]
    genome_pos <- as.matrix(genome_pos[,2, drop = FALSE])

    pmap <- lapply(X = gmap, FUN = function(x) genome_pos[names(x),])
    attr(pmap, "is_x_chr") <- setNames(rep(FALSE, length(pmap)), names(pmap))


  } else {
    prob_list <- lapply(X = data$maps, FUN = function(x) x[[type]]$p1p2$hmm.phase[[1]]$haploprob)
    parents <- setNames(c(1, 2), c(data$data$name.p1, data$data$name.p2))
    ped <- data.frame(id = data$data$screened.data$ind.names, parent1 = data$data$name.p1, parent2 = data$data$name.p2, row.names = NULL)

    # Get the genetic and physical maps
    gmap <- lapply(X = data$maps, FUN = function(x) {
      cm <- c(0, cumsum(imf_h(x[[type]]$p1p2$hmm.phase[[1]]$rf)))
      setNames(cm, row.names(x[[type]]$p1p2$hmm.phase[[1]]$p1))
    })
    names(gmap) <- seq_along(gmap)
    attr(gmap, "is_x_chr") <- setNames(rep(FALSE, length(gmap)), names(gmap))

    # Get the physical map
    pos <- data$data$genome.pos
    genome_pos <- data.frame(marker = names(pos), bp = pos, row.names = NULL)

    row.names(genome_pos) <- genome_pos[[1]]
    genome_pos <- as.matrix(genome_pos[,2, drop = FALSE])

    pmap <- lapply(X = gmap, FUN = function(x) genome_pos[names(x),])
    attr(pmap, "is_x_chr") <- setNames(rep(FALSE, length(pmap)), names(pmap))


  }

  # Create a numeric matrix of covariates for the general combining ability effect (i.e. family)
  ped_covar_mat <- matrix(0, nrow = nrow(ped), ncol = length(parents), dimnames = list(ped$id, names(parents)))
  for (par in names(parents)) {
    ped_covar_mat[, par] <- as.numeric(apply(X = ped, MARGIN = 1, FUN = function(cross) par %in% cross))
  }

  # Reorganize the haplotype (allele) probabilities into a list of arrays
  ped_unique <- unique(ped[,c("parent1", "parent2")])
  # Names of parent haplotypes
  parent_hapl_names <- paste0(rep(names(parents), each = ploidy), "_", rep(paste0("hap", seq_len(ploidy)), length(parents)))

  # Iterate over linkage groups
  allele_probs_list <- list()
  for (i in seq_along(prob_list)) {
    chrom <- as.character(i)
    tmp <- as.matrix(prob_list[[i]])
    # Get the marker names for this lg
    if (class(data) == "mappoly2.consensus.map") {
      mars <- row.names(data$consensus.map[[i]]$ph$G)
      colnames(tmp)[-1:-3] <- mars
    } else {
      mars <- row.names(data$maps[[i]][[type]]$p1p2$hmm.phase[[1]]$p1)
      colnames(tmp) <- c("parent", "ind", "homolog", mars)
      # Get id names
      ids <- data$data$screened.data$ind.names
      # Create rownames for tmp
      tmp_rownames <- paste0(ids[tmp[,"ind"]], ".Par", tmp[,"parent"])
      row.names(tmp) <- tmp_rownames
    }

    # Create an array of haplotype probabilities
    # Dimensions should be indiv x genotypes x markers
    #
    Ind <- ped$id
    nInd <- length(Ind)
    nMar <- length(mars)
    nGeno <- length(parent_hapl_names)

    haplo_prob_array <- array(data = as.numeric(0), dim = c(nInd, nGeno, nMar), dimnames = list(Ind, parent_hapl_names, mars))

    # Use columns of tmp to create indices to place the probabilities in the array
    x_idx <- tmp[,"ind"]
    y_idx <- ((tmp[,"parent"] - 1) * ploidy) + tmp[,"homolog"]
    # Now add the probabilities to the array
    for (j in seq_len(nrow(tmp))) {
      haplo_prob_array[x_idx[j], y_idx[j], ] <- tmp[j,-1:-3]
    }

    # Divide by the ploidy to get the allele fractions
    haplo_prob_array <- haplo_prob_array / ploidy
    # haplo_prob_array <- haplo_prob_array


    # Add the array to the list
    allele_probs_list[[chrom]] <- haplo_prob_array

  }

  # Add attributes
  class(allele_probs_list) <- c("calc_genoprob", "list")
  attr(allele_probs_list, "crosstype") <- ""
  attr(allele_probs_list, "is_x_chr") <- setNames(rep(FALSE, length(allele_probs_list)), names(allele_probs_list))
  attr(allele_probs_list, "alleles") <- parent_hapl_names
  attr(allele_probs_list, "alleleprobs") <- TRUE

  # Return the allele probability matrix, the genetic map, and the physical map
  out <- list(probs = allele_probs_list, gmap = gmap, pmap = pmap, fam.covar = ped_covar_mat)
  return(out)

}










# Convert the mappoly data object to mappoly2 format
mappoly_dat_to_mappoly2_csv <- function(x, path, parent.names = NULL) {
  if (inherits(x, "mappoly.data")) {
    x <- list(x)
  }
  if (is.null(parent.names))
    parent.names <- t(apply(matrix(1:(length(x) * 2), ncol = 2,
                                   byrow = TRUE), 1, function(x) paste0("P", x)))
  assertthat::assert_that(is.matrix(parent.names))
  assertthat::assert_that(ncol(parent.names) == 2 & nrow(parent.names) ==
                            length(x))
  assertthat::assert_that(all(sapply(x, function(x) inherits(x, "mappoly.data"))))
  for (i in 1:length(x)) {
    F1 <- x[[i]]$geno.dose
    F1[F1 == x[[i]]$ploidy + 1] <- NA
    r <- is.null(x[[i]]$seq.ref)
    if (is.null(r))
      r <- rep(NA, length(x[[i]]$mrk.names))
    a <- is.null(x[[i]]$seq.alt)
    if (is.null(r))
      a <- rep(NA, length(x[[i]]$mrk.names))
    w <- data.frame(snp_id = x[[i]]$mrk.names, P1 = x[[i]]$dosage.p1,
                    P2 = x[[i]]$dosage.p2, chrom = x[[i]]$chrom[x[[i]]$mrk.names],
                    genome_pos = x[[i]]$genome.pos[x[[i]]$mrk.names],
                    ref = r, alt = a, F1)
    names(w)[c(2:3)] <- parent.names[i, ]
    fn <- paste0(paste0(parent.names[i, ], collapse = "x"), ".csv")
    if (is.null(path)) {
      file_path <- file.path(getwd(), fn)
    } else {
      file_path <- file.path(path, fn)
    }
    write.csv(w, file = file_path, row.names = FALSE)
  }
}



#' Pull the consensus map
#'
#' @param x A \code{mappoly2.consensus.map} object.
#' @param lg The linkage groups to pull (defaults to all).
#' @param type The type of map to return (can be genetic, physical, or both).
#' @param format The format of the output (can be a table or list format for qtl2).
#'
#' @importFrom qtl map2table
#' @importFrom mappoly2 imf_h
#'
#' @export
#'
pull_consensus_map <- function(x, lg, type = c("genetic", "physical", "both"), format = c("table", "qtl2")) {

  stopifnot(inherits(x, "mappoly2.consensus.map"))
  type <- match.arg(type)
  format <- match.arg(format)
  # group <- match.arg(group)
  if (missing(lg)) {
    lg <- "all"
  }

  # # If type == consensus, check if there are consensus maps
  # if (group == "consensus") {
  #   if (!any(names(x) == "consensus.map")) {
  #     stop("Input object does not have a consensus map.")
  #   } else {

  # Match up LGs
  if (lg == "all") {
    lg <- names(x$consensus.map)
  } else {
    if (is.numeric(lg)) {
      lg <- names(x$consensus.map)[lg]
    } else {
      lg <- intersect(lg, names(x$consensus.map))
    }
  }

  # Get the genetic and physical maps
  gmap <- lapply(X = x$consensus.map[lg], FUN = function(x) {
    cm <- c(0, cumsum(imf_h(x$rf)))
    setNames(cm, row.names(x$ph$G))
  })
  names(gmap) <- as.numeric(sub(pattern = "[a-zA-Z]{1,}", replacement = "", x = names(gmap)))

  # Get the physical map
  genome_pos <- lapply(X = x$individual.maps, FUN = function(x) {
    pos <- x$data$genome.pos
    data.frame(marker = names(pos), bp = pos, row.names = NULL)
  })
  genome_pos <- unique(`row.names<-`(do.call(rbind, genome_pos), NULL))

  # Resolve issues where the same marker has different genome positions; take the first position
  genome_pos <- split(genome_pos, genome_pos$marker)
  genome_pos <- lapply(X = genome_pos, FUN = function(x) x[1,, drop = FALSE])
  genome_pos <- `row.names<-`(do.call(rbind, genome_pos), NULL)

  row.names(genome_pos) <- genome_pos[[1]]
  genome_pos <- as.matrix(genome_pos[,2, drop = FALSE])

  pmap <- lapply(X = gmap, FUN = function(x) genome_pos[names(x),])

  # output format
  if (format == "qtl2") {
    attr(gmap, "is_x_chr") <- setNames(rep(FALSE, length(gmap)), names(gmap))
    attr(pmap, "is_x_chr") <- setNames(rep(FALSE, length(pmap)), names(pmap))

    if (type == "genetic") {
      output <- gmap
    } else if (type == "physical") {
      output <- pmap
    } else {
      output <- list(gmap = gmap, pmap = pmap)
    }

  } else {
    gmap_df <- qtl::map2table(gmap)
    gmap_df <- cbind(gmap_df, marker = row.names(gmap_df))
    row.names(gmap_df) <- NULL
    gmap_df <- gmap_df[c("marker", "chr", "pos")]
    names(gmap_df)[3] <- "cM"

    pmap_df <- qtl::map2table(pmap)
    pmap_df <- cbind(pmap_df, marker = row.names(pmap_df))
    row.names(pmap_df) <- NULL
    pmap_df <- pmap_df[c("marker", "chr", "pos")]
    names(pmap_df)[3] <- "bp"

    if (type == "genetic") {
      output <- gmap_df
    } else if (type == "physical") {
      output <- pmap_df
    } else {
      output <- merge(x = gmap_df, y = pmap_df, by = c("marker", "chr"))
    }
  }

  # Return the output
  return(output)

}



