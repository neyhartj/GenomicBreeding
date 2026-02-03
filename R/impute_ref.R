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
