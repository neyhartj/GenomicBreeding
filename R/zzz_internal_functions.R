#' Internal functions
#'
#' @description
#' Internal package documents or those without clear documentation.
#'
impute_linkage <- function(hd.file, ld.file, out.file, hd.geno, ld.geno, ped, n.core = 1) {
  # Error checking
  if (!missing(hd.file)) {
    stopifnot(file.exists(hd.file))
    if (!grepl(".csv", x = hd.file)) stop("Input 'hd.file' must be a .csv file.")
    hd.geno <- read.csv(file = hd.file, check.names = FALSE)
  } else if (missing(hd.geno)) {
    stop("Either 'hd.file' or 'hd.geno' must be provided.")
  } else {
    # Else validate hd.geno
    stopifnot(is.data.frame(hd.geno))
  }
  # Error checking
  if (!missing(ld.file)) {
    stopifnot(file.exists(ld.file))
    if (!grepl(".vcf", x = ld.file)) stop("Input 'ld.file' must be a .vcf file.")
    vcf_in_ld <- read.vcfR(file = in.file, verbose = FALSE)
    ld.geno <- extract.gt(x = vcf_in_ld)
    ld.geno <- GT2DS(GT = ld.geno, n.core = n.core)
  } else if (missing(ld.geno)) {
    stop("Either 'ld.file' or 'ld.geno' must be provided.")
  } else {
    # Validate ld.geno
    stopifnot(is.matrix(ld.geno))

  }

  stopifnot(is.data.frame(ped))

  # Validate IDs
  low.id <- colnames(ld.geno)
  id <- ped$id[ped$pop > 0]
  if (length(setdiff(id, low.id)) > 0)  stop("Low density file is missing individuals in the pedigree file.")

  if ("bp" %in% colnames(hd.geno)) {
    map <- hd.geno[, c(1, 2, 4)]
    gen_map <- hd.geno[, c(1, 2, 3)]
    hd.geno1 <- hd.geno[, -4]
  } else {
    gen_map <- map <- hd.geno[, 1:3]
    hd.geno1 <- hd.geno
  }
  colnames(hd.geno1) <- replace(colnames(hd.geno1), 1:3, c("marker",  "chrom", "pos"))
  colnames(map) <- c("marker", "chrom", "pos")
  cross_map <- gen_map[,c("chrom", "position")]
  row.names(cross_map) <- gen_map$marker
  cross_map_list <- table2map(tab = cross_map)
  cross_map_list <- lapply(cross_map_list, function(x) rbind(x, x))


  # Subset markers present in the hd.geno data; add NA for markers in the hd.geno
  # file but missing in the ld.geno object
  ld.geno1 <- ld.geno[rownames(ld.geno) %in% hd.geno1$marker, ]
  hd_mars_not_in_ld <- setdiff(hd.geno1$marker, row.names(ld.geno1))
  if (length(hd_mars_not_in_ld) > 0) {
    ld.geno.miss <- matrix(NA, nrow = length(hd_mars_not_in_ld), ncol = ncol(ld.geno1),
                           dimnames = list(hd_mars_not_in_ld, colnames(ld.geno1)))
    ld.geno1 <- rbind(ld.geno1, ld.geno.miss)[hd.geno1$marker, , drop = FALSE]
  }

  # Remove individuals in the ld matrix that are present in the hd matrix
  ld.geno1 <- ld.geno1[, setdiff(colnames(ld.geno1), colnames(hd.geno1)), drop = FALSE]
  geno.inds <- union( colnames(ld.geno1), colnames(hd.geno1)[-1:-3])

  # Create a temporary dir
  if (!dir.exists("tmp"))  dir.create("tmp")
  # Read a fake vcf file in
  blank_vcf <- read.vcfR(file = system.file("extdata/vcf_example_out.vcf.gz", package = "onemap"), verbose = FALSE)
  # Adjust the fix information
  blank_vcf@fix <- as.matrix(cbind(map[,c("chrom", "pos", "marker")], REF = "A", ALT = "T", QUAL = NA, FILTER = "PASS", INFO = NA))
  colnames(blank_vcf@fix)[1:3] <- c("CHROM", "POS", "ID")

  # List of populations that were imputed already
  imputed_pops <- numeric()
  # Vector of unique pops
  all_pops <- setdiff(unique(ped$pop), 0)
  # New high-density matrix
  hd.geno2 <- hd.geno1

  # Iterate over the pops
  for (p in all_pops) {
    # Subset the pedigree for those pops
    ids <- intersect(ped[ped$pop == p, "id"], geno.inds)
    ped_p <- subset_pedigree(ped = ped, ids = ids)

    # Traverse down the pedigree
    ped_p_pops <- setdiff(unique(ped_p$pop), 0)
    for (pp in ped_p_pops) {
      # Has this pop already been imputed?
      if (pp %in% imputed_pops) {
        next
      } else {
        # Subset the pop and its parents
        pp_ped <- subset(ped_p, pop == pp)
        pp_ids <- pp_ped$id
        pp_pars <- unique(unlist(pp_ped[,c("mother", "father")]))

        par_geno <- hd.geno2[,pp_pars, drop = FALSE]
        par_geno_split <- matrix(NA, nrow = nrow(par_geno), ncol = 4)
        par_geno_split[,1:2] <- do.call(rbind, strsplit(x = par_geno[,1], split = "\\|"))
        par_geno_split[,3:4] <- do.call(rbind, strsplit(x = par_geno[,2], split = "\\|"))
        id_geno <- ld.geno1[,pp_ids, drop = FALSE]

        # Determine genotype states
        # 1. markers where both parents are homozygous are uninformative for phasing; set as NA
        mar_both_par_hom <- rowSums(par_geno == "0|0" | par_geno == "1|1") == 2
        id_geno[mar_both_par_hom,] <- NA
        # 2. Assign parent states to each genotype
        id_geno1 <- id_geno
        for (j in seq_len(ncol(id_geno1))) {
          x <- id_geno1[,j]
          nNA <- sum(is.na(x))
          for (m in seq_along(x)) {
            if (is.na(x[m])) {
              next
            } else {
              id_geno_m <- x[m]
              par_geno_m <- par_geno_split[m,]
              if (id_geno_m == 0) {
                w <- paste0(which(par_geno_m == id_geno_m), collapse = "")
                ww <- switch(w, "134" = 5, "234" = 6, "124" = 8, "123" = 7, "13" = 1, "14" = 2, "23" = 3, "24" = 4)
                x[m] <- ifelse(is.null(ww), NA, ww)

              } else if (id_geno_m == 1) {
                w <- paste0(which(par_geno_m == id_geno_m), collapse = "")
                ww <- switch(w, "1" = 5, "2" = 6, "3" = 7, "4" = 8, "13" = 10, "14" = 9, "23" = 9, "24" = 10,
                             "123" = 8, "124" = 7, "134" = 6, "234" = 5)
                x[m] <- ifelse(is.null(ww), NA, ww)

              } else {
                w <- paste0(which(par_geno_m == 1), collapse = "")
                ww <- switch(w, "134" = 5, "234" = 6, "123" = 7, "124" = 8, "13" = 1, "14" = 2, "23" = 3, "24" = 4)
                x[m] <- ifelse(is.null(ww), NA, ww)
              }
            }
          }
          id_geno1[,j] <- x
        }

        # Create a QTL object
        cs <- sim.cross(map = cross_map_list, n.ind = ncol(id_geno1), type = "4way")
        # Set geno
        new_geno <- lapply(cs$geno, FUN = function(xx) {
          xx$data <- t(id_geno1[colnames(xx$data),])
          xx
        })
        cs$geno <- new_geno

        # Run phasing
        cs <- calc.genoprob(cross = cs, error.prob = 0.03)
        cs <- argmax.geno(cross = cs, error.prob = 0.03)

        # Convert to phased genotypes
        phased_ind_geno_probs <- pull.genoprob(cross = cs, rotate = TRUE)
        colnames(phased_ind_geno_probs) <- colnames(id_geno1)
        phased_ind_geno <- pull.argmaxgeno(cross = cs, rotate = TRUE)
        colnames(phased_ind_geno) <- colnames(id_geno1)
        phased_ind_geno <- apply(X = phased_ind_geno, MARGIN = 2, FUN = function(x) {
          x1 <- as.list(x)
          x1[x == 1] <- list(c(1, 3))
          x1[x == 2] <- list(c(1, 4))
          x1[x == 3] <- list(c(2, 3))
          x1[x == 4] <- list(c(2, 4))
          x2 <- mapply(seq_along(x1), x1, FUN = function(i, j) par_geno_split[i,j])
          apply(X = x2, MARGIN = 2, FUN = paste0, collapse = "|")
        })
        row.names(phased_ind_geno) <- row.names(id_geno1)

        hd.geno2 <- cbind(hd.geno2, phased_ind_geno)

        # Run basic imputation


        # Add the pp to the imputed list
        imputed_pops <- c(imputed_pops, pp)

      } # End of if statement

    } # End of loop traversing the pedigree

  } # End of population loop

} # End of function











