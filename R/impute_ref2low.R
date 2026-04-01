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
  library_files <- vcf2ag(vcf = vcf_ref, outdir = tmp_dir, sep.by.chrom = TRUE)
  # convert the input VCF to alpha-genes genotype files
  input_files <- vcf2ag(vcf = vcf_in, outdir = tmp_dir, sep.by.chrom = TRUE)

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
