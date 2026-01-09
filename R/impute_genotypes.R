#' Create a reference set of haplotypes from a high-density set
#'
#' @param in.file Input VCF file.
#' @param out.file Output VCF file.
#' @param mar.max.missing Maximum missingness threshold for markers.
#' @param ind.max.missing Maximum missingness threshold for individuals.
#' @param method Imputation and phasing method ("beagle" for Beagle, "api2" for AlphaPlantImpute2)
#' @param keep.ind A vector of individuals to keep in the VCF file.
#' @param keep.mar A vector of markers to keep in the VCF file.
#' @param params A \code{list} of parameters for beagle.
#' @param beagle.path Path to the executable file for beagle.
#' @param n.core The number of cores for processing.
#'
#' @import vcfR
#'
#' @export
#'
impute_ref <- function(in.file, out.file, mar.max.missing, ind.max.missing, method = c("beagle", "api2"),
                       keep.ind, keep.mar, params, beagle.path, n.core = 1) {
  # Error checking
  stopifnot(file.exists(in.file))
  stopifnot(mar.max.missing >= 0 & mar.max.missing <= 1)
  stopifnot(ind.max.missing >= 0 & ind.max.missing <= 1)
  method <- match.arg(method)

  stopifnot(n.core >= 1)
  n.core <- as.integer(n.core)
  if (missing(params)) {
    if (method == "beagle") {
      params <- list(burnin = 20, iterations = 50, nthreads = n.core)
    } else {
      params <- list(hd_threshold = ind.max.missing, n_haplotypes = 150, error = 0.03, n_sample_rounds = 15, maxthreads = n.core)
    }
  } else {
    if (method == "beagle") {
      params <- params[names(params) != "nthreads"]
      params <- c(params, list(nthreads = n.core))
    } else {
      params <- params[names(params) != "maxthreads"]
      params <- c(params, list(maxthreads = n.core))
    }
  }
  if (!missing(keep.ind)) {
    stopifnot(is.character(keep.ind))
  } else {
    keep.ind <- NULL
  }
  if (!missing(keep.mar)) {
    stopifnot(is.character(keep.mar))
  } else {
    keep.mar <- NULL
  }

  # Read in the vcf
  vcf_in <- read.vcfR(file = in.file, verbose = FALSE)
  # Convert to dosage
  gt <- extract.gt(x = vcf_in, element = "GT")
  geno <- gt
  geno[geno == "0/0"] <- 0
  geno[geno %in% c("0/1", "1/0")] <- 1
  geno[geno == "1/1"] <- 2
  geno <- apply(X = geno, MARGIN = 2, FUN = as.numeric)
  row.names(geno) <- row.names(gt)

  # Filter on missing
  ind_idx <- which(colMeans(is.na(geno)) < ind.max.missing)
  marker_idx <- which(rowMeans(is.na(geno)) < 1)
  geno1 <- geno[marker_idx, ind_idx, drop = FALSE]

  ind_idx <- which(colMeans(is.na(geno1)) < 1)
  marker_idx <- which(rowMeans(is.na(geno1)) < mar.max.missing)
  geno2 <- geno1[marker_idx, ind_idx, drop = FALSE]

  # Remove duplicate genotypes
  if (is.null(keep.ind)) {
    idx <- which(!endsWith(x = colnames(geno2), suffix = ".1"))
  } else {
    idx <- which(colnames(geno2) %in% keep.ind)
  }
  if (is.null(keep.mar)) {
    idy <- seq_len(nrow(geno2))
  } else {
    idy <- which(row.names(geno2) %in% keep.mar)
  }

  geno3 <- geno2[idy, idx, drop = FALSE]

  # Filter on variance
  idx <- which(apply(X = geno3, MARGIN = 1, FUN = sd, na.rm = TRUE) > 0)
  geno4 <- geno3[idx, , drop = FALSE]

  # Create a new vcf with these markers and genotypes
  idx_mar <- which(vcf_in@fix[,"ID"] %in% row.names(geno4))
  idx_geno <- c(1, which(colnames(vcf_in@gt) %in% colnames(geno4)))
  new_vcf <- vcf_in[idx_mar, idx_geno]

  cat("\nReference VCF for imputation has", nrow(new_vcf), "markers and", ncol(new_vcf@gt)-1, "individuals.")

  tmp_dir <- "tmp"
  if (!dir.exists(tmp_dir)) dir.create(tmp_dir)

  # Split on method
  if (method == "beagle") {
    cat("\nRunning imputation and phasing using Beagle...")

    # Save this vcf as a temporary file
    tmp_file <- file.path(tmp_dir, "beagle_input.vcf.gz")
    write.vcf(x = new_vcf, file = tmp_file)

    # Dirname of the base file
    out_dir <- dirname(out.file)
    if (!dir.exists(out_dir)) dir.create(path = out_dir)
    # Remove any extension from the output file
    output <- file.path(out_dir, gsub(pattern = "(.*)(\\..*$)", replacement = "\\1", x = basename(out.file)))

    # Run beagle imputation for the high-density set
    exp <- paste0("java -jar ", beagle.path, " gt=", tmp_file,
                  " out=", output, " ", paste(paste0(names(params), "=", params), collapse = " "))
    stdout <- capture.output(system(command = exp, ignore.stdout = TRUE))

  } else {
    cat("\nRunning imputation and phasing using AlphaPlantImpute2")

    # Alpha Plant Impute2
    # Convert the vcf to alpha genes format
    api2_input_files <- vcf2ag(vcf = new_vcf, outdir = tmp_dir, sep.by.chrom = TRUE)

    # Create parameters for AlphaPlantImpute
    api2_params <- paste0("-", names(params), " ", params)

    # Iterate over the input files
    for (input_file in api2_input_files) {
      input_file_path <- file.path(getwd(), input_file)
      output_file_path <- sub("input", "haplotype_library", input_file)
      output_file_path <- file.path(getwd(), sub(".geno", "", output_file_path))

      exp <- paste0("AlphaPlantImpute2 -createlib -genotypes ", input_file_path, " -out ", output_file_path,
                    " ", paste0(api2_params, collapse = " "))
      stdout <- capture.output(system(command = exp, ignore.stdout = TRUE))

    }

    # Split vcf by chromosome
    marks <- new_vcf@fix[,"ID"]
    markers_by_chrom <- split(marks, new_vcf@fix[,"CHROM"])
    # List the output files
    api2_output_files <- list.files(path = tmp_dir, pattern = ".phase", full.names = TRUE)
    output_mat_list <- list()
    # Iterate over haplotype files
    for (output_file in api2_output_files) {
      df <- read.delim(file = output_file, sep = " ", header = FALSE)
      chrom <- sub(pattern = "(.*_)(chr[0-9]{2,})(\\.phase)", replacement = "\\2", x = output_file)
      names(df) <- c("id", markers_by_chrom[[chrom]])
      df_split <- split(df, df$id)
      df1 <- lapply(X = df_split, FUN = function(x) {
        x1 <- as.matrix(x)[,-1]
        paste0(x1[1,], "|", x1[2,])
      })
      mat <- do.call("rbind", df1)
      row.names(mat) <- names(df1)
      mat <- t(mat)
      row.names(mat) <- names(df)[-1]
      output_mat_list[[chrom]] <- mat

    }

    # Merge the matrices
    gt_new <- do.call(rbind, output_mat_list)[marks,]

    # Create a new VCF object
    vcf_phase <- new_vcf
    vcf_phase@gt <- cbind("FORMAT" = "GT", gt_new)

    meta <- c(vcf_phase@meta, "##Haplotype library built using AlphaPlantImpute2", paste0("##", exp))
    vcf_phase@meta <- meta

    # Save the new vcf file
    output_filename <- file.path(dirname(out.file), paste0(sub(pattern = "\\..*$", replacement = "", x = basename(out.file)), ".vcf.gz"))
    write.vcf(x = vcf_phase, file = output_filename)


  }

  # Remove the temporary folder
  unlink("tmp", recursive = TRUE)

}



#' Impute and phase a low-density offspring set using a high-density phased reference set
#'
#' @param in.file Input VCF file.
#' @param out.file Output VCF file.
#' @param ref.file The phased reference VCF file.
#' @param ped.file The pedigree file or a ped object.
#' @param ind.max.missing Maximum missingness threshold for individuals.
#' @param beagle.params A \code{list} of parameters for beagle.
#' @param api2.params A \code{list} of parameters for AlphaPlantImpute2.
#' @param update.ref Logical. Should the reference set be updated with the low-density phased haplotypes?
#' @param n.core The number of cores for processing.
#'
#' @import vcfR
#' @import polyBreedR
#'
#' @export
#'
impute_ref2low <- function(in.file, out.file, ref.file, ped.file, ind.max.missing, beagle.params, api2.params,
                           update.ref = TRUE, n.core = 1) {

  # Error checking
  if (!inherits(in.file, "vcfR")) {
    stopifnot(file.exists(in.file))
    vcf_in <- read.vcfR(file = in.file, verbose = FALSE)
  } else {
    vcf_in <- in.file
  }
  if (!inherits(ref.file, "vcfR")) {
    stopifnot(file.exists(ref.file))
    vcf_ref <- read.vcfR(file = ref.file, verbose = FALSE)
  } else {
    vcf_ref <- ref.file
  }
  if (missing(ped.file)) {
    ped <- NULL
  } else {
    if (is.data.frame(ped.file)) {
      ped <- ped.file
    } else {
      stopifnot(file.exists(ped.file))
      ped <- get_pedigree2(ped.file = ped.file, trim = FALSE)
    }
    stopifnot(c("id", "mother", "father") %in% names(ped))
  }
  stopifnot(ind.max.missing >= 0 & ind.max.missing <= 1)
  stopifnot(is.logical(update.ref))
  stopifnot(n.core >= 1)
  n.core <- as.integer(n.core)
  if (missing(api2.params)) {
    api2.params <- list(hd_threshold = ind.max.missing, n_haplotypes = 150, error = 0.03, n_sample_rounds = 15, maxthreads = n.core)
  } else {
    api2.params <- api2.params[names(api2.params) != "maxthreads"]
    api2.params <- c(api2.params, list(maxthreads = n.core))
  }
  if (missing(beagle.params)) {
    beagle.params <- list(burnin = 20, iterations = 50, nthreads = n.core)
  } else {
    beagle.params <- beagle.params[names(beagle.params) != "nthreads"]
    beagle.params <- c(beagle.params, list(nthreads = n.core))
  }
  # Create parameters for AlphaPlantImpute
  api2_params <- paste0("-", names(api2.params), " ", api2.params)

  ## Trim the input vcf file for individuals (not in the refererence file) and markers (only those in the reference file)
  ref_inds <- colnames(vcf_ref@gt)[-1]
  ref_marks <- vcf_ref@fix[,"ID"]
  idx_mar <- which(vcf_in@fix[,"ID"] %in% ref_marks)
  idx_ind <- which(! colnames(vcf_in@gt) %in% ref_inds)
  vcf_in <- vcf_in[idx_mar, c(1, idx_ind)]

  # Subset individuals based on missing data
  gt <- extract.gt(x = vcf_in, element = "GT")
  gt[gt == "./."] <- as.character(NA)
  idx <- which(colMeans(is.na(gt)) <= ind.max.missing)
  vcf_in <- vcf_in[,idx]
  vcf_inds <- colnames(vcf_in@gt)[-1]

  cat("Low-density VCF for imputation has", nrow(vcf_in), "markers and", ncol(vcf_in@gt)-1, "individuals.")

  if (!is.null(ped)) {
    # Create a pedigree founder file
    ped1 <- subset_pedigree(ped = ped, ids = vcf_inds)
    # Make sure all ids and parents are in the ref vcf or input vcf
    ped_ids <- na.omit(unique(unlist(ped1[,c("id", "mother", "father")])))
    if (!all(ped_ids %in% union(ref_inds, vcf_inds))) stop("Pedigree individuals and founders are not in the input or reference VCF file.")
    ped1 <- ped1[,c("id", "mother", "father")]
    # Save it
    write.table(x = ped1, file = file.path(tmp_dir, "alphaplantimpute_founders.txt"), quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
    ped_arg <- paste0("-founder ", file.path(tmp_dir, "alphaplantimpute_founders.txt"))

  } else {
    ped_arg <- ""

  }

  # Create temporary directory
  tmp_dir <- "tmp"
  if (!dir.exists(tmp_dir)) dir.create(tmp_dir)

  # Convert the reference file to phased library file
  library_files <- berryBreedR:::vcf2ag(vcf = vcf_ref, outdir = tmp_dir, sep.by.chrom = TRUE)
  # convert the input VCF to alpha-genes genotype files
  input_files <- berryBreedR:::vcf2ag(vcf = vcf_in, outdir = tmp_dir, sep.by.chrom = TRUE)

  # Iterate over chromosomes and run imputation
  chrnames <- unique(vcf_ref@fix[,"CHROM"])
  for (chrom in chrnames) {

    libfile <- grep(pattern = chrom, x = library_files, value = TRUE)
    genofile <- grep(pattern = chrom, x = input_files, value = TRUE)
    output_file_path <- sub("input", "output", genofile)
    output_file_path <- file.path(getwd(), sub(".geno", "", output_file_path))

    exp <- paste0("AlphaPlantImpute2 -impute -libphase ", libfile, " -genotypes ", genofile, " -out ", output_file_path,
                  " ", paste0(api2_params, collapse = " "), " ", ped_arg)
    stdout <- capture.output(system(command = exp, ignore.stdout = TRUE))

  }

  ## Convert the imputed genotypes to a VCF file
  # Split vcf by chromosome
  marks <- vcf_in@fix[,"ID"]
  markers_by_chrom <- split(marks, vcf_in@fix[,"CHROM"])
  # List the output files
  api2_output_files <- list.files(path = tmp_dir, pattern = ".genotypes", full.names = TRUE)
  output_mat_list <- list()
  # Iterate over haplotype files
  for (output_file in api2_output_files) {
    df <- read.delim(file = output_file, sep = " ", header = FALSE)
    chrom <- sub(pattern = "(alphaplantimpute_output_)(.*chr[0-9]{2,})(\\.genotypes)", replacement = "\\2", x = basename(output_file))
    colnames(df) <- c("id", markers_by_chrom[[chrom]])
    mat <- t(df[-1])
    colnames(mat) <- df$id
    mat[mat == 0] <- "0/0"
    mat[mat == 1] <- "0/1"
    mat[mat == 2] <- "1/1"
    output_mat_list[[chrom]] <- mat
  }

  # Merge the matrices
  gt_new <- do.call(rbind, output_mat_list)[marks,]

  # Create a new VCF object
  vcf_imputed <- vcf_in
  vcf_imputed@gt <- cbind("FORMAT" = "GT", gt_new)

  # Save the imputed genotypes
  out.file1 <- paste0(out.file, ".vcf.gz")
  write.vcf(x = vcf_imputed, file = out.file1)
  # Notify
  cat("\nImputed genotypes written to:", out.file1)

  # Update the phased haplotypes if called
  if (update.ref) {
    vcf_new_ref <- vcf_ref
    vcf_new_ref@gt <- cbind(FORMAT = "GT", vcf_ref@gt[,-1], vcf_imputed@gt[,-1])
    # Save the VCF as a new file
    beagle_input_file <- file.path(tmp_dir, "beagle_haplotype_update_genotypes.vcf.gz")
    write.vcf(x = vcf_new_ref, file = beagle_input_file)

    outout_updated_libfile <- sub(pattern = "library", replacement = "library_updated", x = ref.file)
    outout_updated_libfile <- sub(pattern = "(.*)(\\.vcf)(.*)", replacement = "\\1", x = outout_updated_libfile)
    # Run beagle imputation for the high-density set
    exp <- paste0("java -jar /opt/beagle/beagle.06Aug24.a91.jar gt=", beagle_input_file,
                  " out=", outout_updated_libfile, " ", paste(paste0(names(beagle.params), "=", beagle.params), collapse = " "))
    stdout <- capture.output(system(command = exp, ignore.stdout = TRUE))
  }

  # Remove the temporary folder
  unlink("tmp", recursive = TRUE)

  # Notify
  cat("\nUpdated reference genotypes written to:", paste0(outout_updated_libfile, ".vcf.gz"))

}


# Function to impute from high-density to low-density using linkage analysis
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




#' Impute markers using rrBLUP
#'
#' @description
#' Impute markers using rrBLUP. This function is borrowed from the \code{polyBreedR} package.
#'
#'
#' @param geno A data.frame of genotype calls, where the first three columns are "marker," "chrom,", and "position."
#' Subsequent columns are genotype calls for individuals.
#' @param ploidy Integer ploidy number.
#' @param method Imputation method: the population mean ("pop"), expectation-maximization ("EM"), or random forest ("RF").
#' @param ind.max.missing Maximum missingness threshold for individuals.
#' @param snp.max.missing Maximum missingness threshold for markers.
#' @param params A \code{list} of parameters.
#' @param n.core The number of cores for processing.
#'
#' @import polyBreedR
#' @import rrBLUP
#'
#' @export
#'
impute_geno <- function(geno, ploidy, method = c("pop", "EM", "RF"), ind.max.missing, snp.max.missing,
                        params = NULL, n.core = 1) {

  stopifnot(is.data.frame(geno))
  cols.match <- names(geno)[1:3] == c("marker", "chrom", "position")
  if (any(!cols.match)) {
    stop("The first 3 columns of 'geno' should be 'marker', 'chrom', and 'position'.")
  }
  method <- match.arg(method)
  if (!is.null(params) && !is.list(params))
    stop("Use a list for argument params")
  if (method == "EM")
    params1 <- list(tol = 0.02)
  if (method == "RF")
    params1 <- list(ntree = 100, nflank = 100, tol = 0.02)
  if (!is.null(params)) {
    if (method == "EM")
      stopifnot(names(params) == "tol")
    if (method == "RF")
      stopifnot(names(params) %in% c("ntree", "nflank",
                                     "tol"))
    params1[names(params)] <- params
  }
  impute.mode <- function(x) {
    miss <- which(is.na(x))
    if (length(miss) > 0) {
      x[miss] <- as.integer(names(which.max(table(x))))
    }
    return(x)
  }
  impute.mean <- function(x) {
    miss <- which(is.na(x))
    if (length(miss) > 0) {
      x[miss] <- mean(x[-miss])
    }
    return(x)
  }
  impute.RF <- function(y, x, ntree) {
    miss <- which(is.na(y))
    n.miss <- length(miss)
    if (n.miss > 0) {
      ans <- suppressWarnings(randomForest(x = x[-miss,
      ], y = y[-miss], xtest = x[miss, , drop = FALSE],
      ntree = ntree))
      y[miss] <- ans$test$predicted
    }
    if (is.factor(y)) {
      return(as.integer(as.character(y)))
    }
    else {
      return(y)
    }
  }
  geno_orig <- geno
  geno <- as.matrix(geno_orig[,-1:-3])
  row.names(geno) <- all.marks <- geno_orig$marker
  all.marks <- row.names(geno)
  n <- ncol(geno)
  p <- nrow(geno)
  ix <- which(apply(geno, 1, sd, na.rm = T) > 0)
  nd <- nrow(geno) - length(ix)
  if (nd > 0)  cat(sub("X", nd, "Removed X markers without genetic variance\n"))
  iu <- which(apply(geno[ix, ], 1, function(x) { sum(is.na(x)) })/n <= snp.max.missing)
  nd2 <- length(ix) - length(iu)
  if (nd2 > 0)  cat(sub("X", nd2, "Removed X markers due to missing data\n"))
  iv <- which(apply(geno[ix, ], 2, function(x) { sum(is.na(x)) })/p <= ind.max.missing)
  ni <- n - length(iv)
  if (nd2 > 0)  cat(sub("X", ni, "Removed X individuals due to missing data\n"))


  geno <- geno[ix[iu], iv, drop = FALSE]
  m <- nrow(geno)
  marks <- rownames(geno)
  ik <- match(marks, all.marks)
  if (method == "pop") {
    geno.imp <- apply(geno, 1, impute.mode)
  }
  if (method %in% c("EM", "RF")) {
    ans <- A.mat(t(geno)/(ploidy/2) - 1, impute.method = "EM",
                 n.core = n.core, return.imputed = TRUE, min.MAF = 0,
                 tol = params1$tol)
    geno.imp <- apply(ans$imputed, 2, function(x) {
      ifelse(abs(x) > 1, 1 * sign(x), x)
    })
    digits <- 0
    geno.imp <- round((geno.imp + 1) * ploidy/2, digits)
  }
  if (method == "RF") {
    geno.imp.old <- geno.imp
    if (n.core > 1) {
      cl <- makeCluster(n.core)
      clusterExport(cl = cl, varlist = NULL)
      geno.imp <- parApply(cl, array(1:m), MARGIN = 1,
                           function(k) {
                             y <- geno[k, ]
                             if (geno.key == "GT")
                               y <- factor(y)
                             impute.RF(y = y, x = geno.imp.old[, setdiff(c(max(1, k - params1$nflank):min(m, k + params1$nflank)),k), drop = FALSE], ntree = params1$ntree)
                           })
      stopCluster(cl)
    } else {
      geno.imp <- apply(array(1:m), MARGIN = 1, function(k) {
        y <- geno[k, ]
        if (geno.key == "GT")
          y <- factor(y)
        impute.RF(y = y, x = geno.imp.old[, setdiff(c(max(1,
                                                          k - params1$nflank):min(m, k + params1$nflank)),
                                                    k), drop = FALSE], ntree = params1$ntree)
      })
    }
  }

  geno.imp <- t(geno.imp)
  row.names(geno.imp) <- marks

  # Merge with snp info
  snp.info <- geno_orig[ik, c("marker", "chrom", "position")]
  geno.imp.df <- cbind(snp.info, geno.imp)

  # Return
  return(geno.imp.df)
}






