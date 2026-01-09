#' Merge VCF files together
#'
#' @param vcf.list A list of \code{vcfR} objects.
#'
#' @return A \code{vcfR} object with merged variant data.
#'
#' @import vcfR
#'
#' @export
#'
#'
merge_vcfs <- function(vcf.list) {

  stopifnot(all(sapply(vcf.list, inherits, "vcfR")))

  complements = c("A" = "T", "C" = "G", "T" = "A", "G" = "C")

  # Get the marker names in each VCF
  marker_names_list <- lapply(X = vcf.list, FUN = function(x) x@fix[,"ID"])
  # Get unique
  marker_names <- unique(unlist(marker_names_list))
  # Get all sample names
  sample_names_list <- lapply(X = vcf.list, FUN = function(x) colnames(x@gt)[-1])
  # Index of where the sample names will be in the larger GT matrix
  sample_idy_list <- list()
  for (i in seq_along(sample_names_list)) {
    if (i == 1) {
      sample_idy_list[[i]] <- seq_along(sample_names_list[[i]])
    } else {
      sample_idy_list[[i]] <- max(sample_idy_list[[i-1]]) + seq_along(sample_names_list[[i]])
    }
  }
  sample_names <- unlist(sample_names_list)

  # Create an empty GT matrix
  gt <- matrix(as.character(NA), nrow = length(marker_names), ncol = length(sample_names),
               dimnames = list(marker_names, sample_names))

  # Iterate over each vcf
  for (i in seq_along(vcf.list)) {
    vcf <- vcf.list[[i]]
    # get the GT
    gt_i <- vcf@gt[,-1]
    mar_names <- vcf@fix[,"ID"]
    # Add it to the larger GT
    idy <- sample_idy_list[[i]]
    gt[mar_names, idy] <- gt_i
  }

  # Replace NAs
  gt[is.na(gt)] <- "./.:0,0:0"

  # Get the FIX information for each vcf; merge
  fix_list <- lapply(X = vcf.list, FUN = getFIX, getINFO = TRUE)
  fix_all <- NULL
  for (i in seq_along(fix_list)) {
    fix <- fix_list[i]
    if (i == 1) {
      fix_all <- as.data.frame(fix)[,c("CHROM", "POS", "ID", "REF", "ALT")]
      names(fix_all)[4:5] <- c("REF.1", "ALT.1")
    } else {
      fix_add <- as.data.frame(fix)[,c("CHROM", "POS", "ID", "REF", "ALT")]
      names(fix_add)[4:5] <- paste0(c("REF", "ALT"), ".", i)
      fix_all <- merge(x = fix_all, y = fix_add, by = c("CHROM", "POS", "ID"), all = TRUE)
    }
  }
  # Sort
  fix_all <- fix_all[order(fix_all$ID), ]

  # Find SNPs where the reference and alternate alleles do not match
  ref_no_match <- apply(X = fix_all[,grep(pattern = "REF", x = colnames(fix_all))], MARGIN = 1, FUN = function(x) length(unique(na.omit(x))) > 1)
  alt_no_match <- apply(X = fix_all[,grep(pattern = "ALT", x = colnames(fix_all))], MARGIN = 1, FUN = function(x) length(unique(na.omit(x))) > 1)

  # Iterate over ref no match
  ref_no_match_idx <- which(ref_no_match)
  for (ii in ref_no_match_idx) {
    refs <- unlist(fix_all[ii, grep(pattern = "REF", x = colnames(fix_all))])
    refs <- refs[!is.na(refs)]
    ref_most_freq <- names(which.max(table(refs)))
    # Are all of the refs either ref most freq or its complement?
    if (all(refs %in% c(ref_most_freq, complements[ref_most_freq]))) {
      # Change all refs
      fix_all[ii, grep(pattern = "REF", x = colnames(fix_all))] <- ref_most_freq
      # Remove ii from the ref_no_match
      ref_no_match_idx <- setdiff(ref_no_match_idx, ii)
    } # Otherwise remove the SNP
  }

  # Do the same thing for alt
  # Iterate over ref no match
  alt_no_match_idx <- which(alt_no_match)
  for (ii in alt_no_match_idx) {
    alts <- unlist(fix_all[ii, grep(pattern = "ALT", x = colnames(fix_all))])
    alts <- alts[!is.na(alts)]
    alt_most_freq <- names(which.max(table(alts)))
    # Are all of the alts either alt most freq or its complement?
    if (all(alts %in% c(alt_most_freq, complements[alt_most_freq]))) {
      # Change all alts
      fix_all[ii, grep(pattern = "ALT", x = colnames(fix_all))] <- alt_most_freq
      # Remove ii from the ref_no_match
      alt_no_match_idx <- setdiff(alt_no_match_idx, ii)
    } # Otherwise remove the SNP
  }

  # Index of SNPs to remove
  allele_no_match <- union(ref_no_match_idx, alt_no_match_idx)

  cat("\nNumber of SNPs removed due to inconsistent REF/ALT alleles:", length(allele_no_match))

  fix_all1 <- fix_all[setdiff(seq_len(nrow(fix_all)), allele_no_match), ]
  # Fix for the consensus ref/alt
  fix_all2 <- fix_all1[,c("CHROM", "POS", "ID")]
  fix_all2$REF <- apply(X = fix_all1[, grepl(pattern = "REF", x = colnames(fix_all1))], MARGIN = 1, FUN = function(x) unique(na.omit(x)))
  fix_all2$ALT <- apply(X = fix_all1[, grepl(pattern = "ALT", x = colnames(fix_all1))], MARGIN = 1, FUN = function(x) unique(na.omit(x)))
  fix_all2$QUAL <- as.character(NA)
  fix_all2$FILTER <- as.character(NA)

  gt1 <- gt[fix_all2$ID, , drop = FALSE]
  colnames(gt1) <- make.unique(colnames(gt1))

  # Recalculate summary statistics
  AD <- apply(X = gt1, MARGIN = 1, FUN = function(x) sapply(strsplit(x = x, split = ":"), "[[", 2))
  AD <- t(AD)
  ref <- masplit(myMat = AD, delim = ",", count = 0, record = 1, sort = 0)
  alt <- masplit(myMat = AD, delim = ",", count = 0, record = 2, sort = 0)
  DP <- ref + alt
  NS <- rowSums(DP > 0)
  AF <- round(rowSums(alt) / rowSums(DP), 3)
  AF[is.na(AF)] <- "."
  DP.AVG <- round(rowMeans(DP), 1)

  info <- paste0("DP=", rowSums(DP), ";NS=", NS, ";DP.AVG=", DP.AVG, ";AF=", AF)

  # Add info flag for target/off-target SNPs
  fix_target_list <- lapply(fix_list, function(x) {
    x1 <- x[,c("ID", "INFO")]
    info_split <- strsplit(x = x[,"INFO"], split = ";")
    target <- sapply(X = info_split, grep, pattern = "TYPE", value = TRUE)
    cbind(ID = x1[,"ID"], INFO = target)
  })
  target_mat <- do.call(rbind, fix_target_list)
  target_mat <- unique(target_mat)
  # In the case of duplicates, choose "target"
  target_mat <- target_mat[order(target_mat[,"INFO"], decreasing = TRUE), ]
  target_mat <- target_mat[!duplicated(target_mat[, "ID"]), ]
  target_mat <- target_mat[order(target_mat[,"ID"]), ]

  idx <- match(x = fix_all2[,"ID"], table = target_mat[,"ID"])
  info <- paste0(info, ";", target_mat[idx, "INFO"])

  fix_all2$INFO <- info

  gt1 <- cbind(FORMAT = "GT:AD:DP", gt1)

  # Assemble the VCF
  vcf_out <- vcf.list[[1]]
  # Sort again
  fix_all2 <- fix_all2[order(fix_all2$ID), ]
  gt1 <- gt1[fix_all2$ID, ]

  vcf_out@fix <- as.matrix(fix_all2)
  vcf_out@gt <- gt1

  # Return!
  return(vcf_out)

}






#' Merge VCF files using allele depth information only
#'
#' @param file.list A character vector of VCF file paths
#' @param run.ids Vector of character defining the run ID for each file
#' @param n.core Number of processors for parallel computation.
#'
#' @import vcfR
#' @import parallel
#'
#' @export
#'
merge_count_vcfs <- function(file.list, run.ids = NULL, n.core = 1) {

  if (!all(endsWith(x = file.list, ".vcf") | endsWith(x = file.list, ".vcf.gz"))) stop("Not all files in 'file.list' are VCF files.")

  # Read in the VCF files
  cat("Reading in the VCF files...\n")
  vcf_list <- lapply(X = file.list, FUN = read.vcfR, verbose = FALSE)

  cat(length(vcf_list), "VCF files passed with...\n")
  cat(paste0("[", paste0(sapply(vcf_list, nrow), collapse = ", "), "]"), "variants and\n")
  cat(paste0("[", paste0(sapply(vcf_list, function(x) ncol(x@gt)-1), collapse = ", "), "]"), "samples.")

  # Get the FIX information for each vcf; merge
  fix_list <- lapply(X = vcf_list, FUN = getFIX, getINFO = TRUE)
  fix_all <- NULL
  for (i in seq_along(fix_list)) {
    fix <- fix_list[i]
    if (i == 1) {
      fix_all <- as.data.frame(fix)[,c("CHROM", "POS", "ID", "REF", "ALT")]
    } else {
      fix_all <- merge(x = fix_all, y = as.data.frame(fix)[,c("CHROM", "POS", "ID", "REF", "ALT")],
                       by = c("CHROM", "POS", "ID"), all = TRUE)
    }
  }
  # Sort
  fix_all <- fix_all[order(fix_all$ID), ]

  # Find SNPs where the reference alleles do not match
  ref_no_match <- apply(X = fix_all[,grep(pattern = "REF", x = colnames(fix_all))], MARGIN = 1, FUN = function(x) length(unique(na.omit(x))) > 1)
  fix_all_ref_no_match <- fix_all[ref_no_match, ]
  cat("\n\nNumber of SNPs removed due to inconsistent REF alleles:", sum(ref_no_match))

  fix_all1 <- fix_all[!ref_no_match, ]

  cat("\n\nMerging read count information in the VCF files...")


  ## Merge AD information
  ad_list <- lapply(X = vcf_list, FUN = extract.gt, element = "AD", IDtoRowNames = TRUE)
  # Get reference allele counts
  ad_ref_list <- lapply(ad_list, masplit, delim = ",", record = 1, sort = 0)
  # Get the first alternate allele count
  ad_alt_list <- lapply(ad_list, masplit, delim = ",", record = 2, sort = 0)

  # Create a full reference allele matrix
  full_ref_ad_mat <- matrix(0, nrow = nrow(fix_all1), ncol = sum(sapply(ad_ref_list, ncol)),
                            dimnames = list(fix_all1$ID, unlist(sapply(ad_ref_list, colnames))))
  # Fill in this matrix
  for (ad in ad_ref_list) {
    ad1 <- ad[intersect(row.names(full_ref_ad_mat), row.names(ad)), , drop = FALSE]
    idx <- match(x = row.names(ad1), table = row.names(full_ref_ad_mat), nomatch = 0)
    idy <- match(x = colnames(ad1), table = colnames(full_ref_ad_mat), nomatch = 0)
    full_ref_ad_mat[idx, idy] <- ad1
  }

  # Do the same thing for the "first" alternate allele
  # Create a full reference allele matrix
  full_alt1_ad_mat <- matrix(0, nrow = nrow(fix_all1), ncol = sum(sapply(ad_alt_list, ncol)),
                             dimnames = list(fix_all1$ID, unlist(sapply(ad_alt_list, colnames))))
  # Fill in this matrix
  for (ad in ad_alt_list) {
    ad1 <- ad[intersect(row.names(full_alt1_ad_mat), row.names(ad)), , drop = FALSE]
    idx <- match(x = row.names(ad1), table = row.names(full_alt1_ad_mat), nomatch = 0)
    idy <- match(x = colnames(ad1), table = colnames(full_alt1_ad_mat), nomatch = 0)
    full_alt1_ad_mat[idx, idy] <- ad1
  }



  # The end result will be an m-length list of n x a of allele depth information at each site
  # Create the empty list
  site_ad_list <- vector("list", length = nrow(fix_all1))
  # Which sites are multiallelic?
  fix_all_alleles <- as.matrix(fix_all1[,-1:-3])
  fix_ref_alleles <- fix_all_alleles[,grepl(pattern = "REF", x = colnames(fix_all_alleles)), drop = FALSE]
  fix_alt_alleles <- fix_all_alleles[,grepl(pattern = "ALT", x = colnames(fix_all_alleles)), drop = FALSE]
  is_multiallelic <- apply(X = fix_all_alleles, MARGIN = 1, FUN = function(x) length(unique(na.omit(unlist(strsplit(x = x, split = ","))))) > 2)

  names(site_ad_list) <- fix_all1$ID

  # Iterate over sites
  for (i in seq_along(site_ad_list)) {
    id <- names(site_ad_list)[i]
    # Get the reference allele
    ref_allele <- unique(na.omit(fix_ref_alleles[i,]))

    # Get the reference allele count and create a matrix
    ref_allele_count <- list(full_ref_ad_mat[id, , drop = FALSE])
    names(ref_allele_count) <- ref_allele

    # Get the alternate allele counts
    # Is the site multi-allelic?
    if (is_multiallelic[i]) {
      alt_allele <- fix_alt_alleles[i,]
      alt_allele_list <- strsplit(x = alt_allele, ",")
      alt_alleles <- unique(na.omit(unlist(alt_allele_list)))
      # How many alternate alleles?
      nAlt <- length(alt_alleles)
      # Create a list of matrices for those alleles
      alt_allele_count <- replicate(n = nAlt,
                                    expr = matrix(0, nrow = 1, ncol = ncol(full_alt1_ad_mat), dimnames = list(id, colnames(full_alt1_ad_mat))),
                                    simplify = FALSE)
      names(alt_allele_count) <- alt_alleles

      # Iterate over the alt_allele_list
      for (j in seq_along(alt_allele_list)) {
        alts <- alt_allele_list[[j]]
        # Get the ad matrix for j
        ad_j <- ad_list[[j]]
        # Skip if id is not in ad_j
        if (id %in% row.names(ad_j)) {
          # Iterate over alternate alleles
          for (k in seq_along(alts)) {
            alt <- alts[k]
            ad_k <- masplit(myMat = ad_j[id, , drop = FALSE], delim = ",", sort = 0, record = k + 1)
            alt_allele_count[[alt]][,colnames(ad_j)] <- ad_k
          }
        }
      }

    } else {
      alt_allele <- unique(na.omit(fix_alt_alleles[i,]))
      # Simply use the alt1 ad matrix
      alt_allele_count <- list(full_alt1_ad_mat[id, , drop = FALSE])
      names(alt_allele_count) <- alt_allele

    }

    # Create the array and add to the list
    site_ad_list[[i]] <- do.call(abind::abind, list(c(ref_allele_count, alt_allele_count), along = 3))

  } # Close the ID loop

  ## Construct the VCF
  #
  # First build the GT
  new_ad_mat <- lapply(X = site_ad_list, FUN = function(x) x[,,])
  new_dp <- lapply(X = new_ad_mat, rowSums)
  new_dp <- do.call(rbind, new_dp)
  new_ad <- lapply(X = new_ad_mat, FUN = function(x) apply(x, 1, paste0, collapse = ","))
  new_ad <- do.call(rbind, new_ad)
  new_ra <- lapply(X = new_ad_mat, function(x) x[,1])
  new_ra <- do.call(rbind, new_ra)
  new_gt <- matrix(paste0(new_dp, ":", new_ra, ":", new_ad), nrow = nrow(new_dp), ncol = ncol(new_dp))
  dimnames(new_gt) <- dimnames(new_dp)
  new_gt <- cbind(FORMAT = "DP:RA:AD", new_gt)

  # Make sample names unique
  colnames(new_gt) <- make.unique(colnames(new_gt))


  # Next build the FIX
  new_fix <- as.matrix(fix_all1[,1:3])
  allele_names <- lapply(X = site_ad_list, FUN = function(x) dimnames(x)[[3]])
  alts <- sapply(X = allele_names, FUN = "[", -1)
  new_fix <- cbind(new_fix, REF = sapply(X = allele_names, FUN = "[[", 1), ALT = sapply(alts, paste0, collapse = ","),
                   QUAL = as.character(NA), FILTER = as.character(NA))
  # Calculate terms for DP, ADS, and TYPE
  INFO_dp <- rowSums(new_dp)
  INFO_ads <- paste0(rowSums(new_ra), ",", INFO_dp - rowSums(new_ra))
  INFO_type <- lapply(X = fix_list, FUN = function(x) x[,c("ID", "INFO")])
  INFO_type <- lapply(X = INFO_type, FUN = function(x) cbind(ID = x[,"ID"], TYPE = ifelse(grepl(pattern = "offtarget", x = x[,"INFO"]), "TYPE=offtarget", "TYPE=target")))
  INFO_type <- do.call(rbind, INFO_type)
  INFO_type <- unique(INFO_type)
  INFO_type <- INFO_type[match(x = new_fix[,"ID"], table = INFO_type[,"ID"]), "TYPE"]
  new_fix <- cbind(new_fix, INFO = paste0("DP=", INFO_dp, ";ADS=", INFO_ads, ";", INFO_type))
  row.names(new_fix) <- NULL

  vcf_out <- vcf_list[[1]]
  vcf_out@fix <- new_fix
  vcf_out@gt <- new_gt

  cat("done!")
  # Return the VCF
  return(vcf_out)

}


# A function to call SNPs using updog from a VCF dosage file
#
# vcf.in - a vcfR object
# out.file - path to the filename for the output VCF
# min.ad.sample - minimum AD of non-zero samples to count towards keeping a site
# min.sample.site - minumum number of samples with >= min.ad.sample to keep a site
# keep.all.target.snps - If TRUE, target SNPs are retained regardless of ad and samples; if FALSE,
# target SNPs will be filtered like any other SNP
# n.core - number of computing cores
#
geno_call_counts <- function(vcf.in, out.file, min.ad.sample = 0, min.sample.site = 0, keep.all.target.snps = TRUE, n.core = 1) {
  require(updog)
  require(vcfR)
  require(parallel)
  require(pbapply)

  stopifnot(inherits(vcf.in, "vcfR"))
  stopifnot(is.character(out.file))
  stopifnot(min.ad.sample >= 0)
  stopifnot(min.sample.site >= 0)
  stopifnot(is.logical(keep.all.target.snps))
  stopifnot(n.core >= 1)

  # Get fixed information
  fixed <- getFIX(x = vcf.in, getINFO = TRUE)
  # If keep.all.target.snps is TRUE, INFO must have a type
  if (keep.all.target.snps) {
    info <- getINFO(x = vcf.in)
    which_target <- grep(pattern = "TYPE=target", x = info, ignore.case = FALSE)
    which_offtarget <- grep(pattern = "TYPE=offtarget", x = info, ignore.case = FALSE)
    if (length(which_target) == 0 & length(which_offtarget) == 0) stop("No target or off-target SNPs. Please make sure that TYPE=[target,offtarget] is specified in the INFO field.")

  } else {
    which_target <- numeric()
    which_offtarget <- seq_len(nrow(fixed))

  }

  # Get allele depth information
  ad <- extract.gt(x = vcf.in, element = "AD")
  # Get the reference matrix
  refmat <- masplit(myMat = ad, delim = ",", record = 1, sort = 0)

  # Get the alternate matrix
  altmat <- refmat
  altmat[] <- 0

  # Get biallelic markers
  idx <- which(is.biallelic(x = vcf.in))
  altmat_biallelic <- masplit(myMat = ad[idx,,drop = FALSE], record = 2, sort = 0)
  altmat[idx, ] <- altmat_biallelic

  # Among multi-allelic sites, which SNPs pass the QC
  idx <- which(!is.biallelic(x = vcf.in))
  # Iterate
  for (i in seq_along(idx)) {
    ii <- idx[i]
    n_alt <- sapply(strsplit(fixed[ii,"ALT"], ","), length)
    ad_alt_pass <- logical(n_alt)
    for (j in seq_len(n_alt)) {
      ad_j <- masplit(myMat = ad[ii,,drop = FALSE], record = j+1, sort = 0)
      ad_alt_pass[j] <- sum(ad_j >= min.ad.sample) >= min.sample.site
    }

    # If only one alternate allele passes the test, add it the altmat
    if (sum(ad_alt_pass) == 1) {
      altmat[ii,] <- masplit(myMat = ad[ii,,drop = FALSE], record = seq_len(n_alt)[which(ad_alt_pass)]+1, sort = 0)
    } else {
      altmat[ii,] <- as.numeric(NA) # These will be removed from the refmat too
    }
  }

  # Find altmat rows that are NA; keep targets
  which_keep <- sort(union(which_target, which(rowMeans(is.na(altmat)) != 1)))
  mar_keep <- row.names(refmat)[which_keep]
  # Tell the user how many variants were removed due to multiallelism
  n_rm <- length(setdiff(seq_len(nrow(altmat)), which_keep))
  cat("\nNumer of multi-allelic SNPs removed:", n_rm)

  # Now determine which sites to keep based on allele and sample depth
  refmat <- refmat[which_keep, ]
  altmat <- altmat[which_keep, ]
  sizemat <- refmat + altmat

  which_keep2 <- which(rowSums(sizemat > min.ad.sample) > min.sample.site)
  mar_keep2 <- row.names(sizemat)[which_keep2]
  mar_keep2 <- sort(union(row.names(refmat)[which_target], mar_keep2))

  n_rm <- length(mar_keep) - length(mar_keep2)
  cat("\nNumer of low-depth and low-sample SNPs removed:", n_rm)

  refmat <- refmat[mar_keep2, ]
  sizemat <- sizemat[mar_keep2, ]

  # Parallelize flexdog
  # Create a function over which to parallelize a list
  run_flexdog <- function(mat.list, ploidy = 2, model = "norm") {
    iidx <- seq_len(nrow(mat.list$rmat))
    snp_names <- row.names(mat.list$rmat)
    ind_names <- colnames(mat.list$rmat)
    snpdf <- NULL
    inddf <- NULL
    # Iterate
    for (i in iidx) {
      out <- updog::flexdog(refvec = mat.list$rmat[i,], sizevec = mat.list$smat[i,], ploidy = ploidy, model = model, verbose = FALSE)
      # Add to the lists
      snpdf_i <- c(snp = snp_names[i],
                   out[c("bias", "seq", "od", "prop_mis", "num_iter", "llike")],
                   ploidy = ploidy, model = model,
                   setNames(replicate(4, NA, simplify = FALSE), c("p1ref", "p1size", "p2ref", "p2size")),
                   setNames(as.list(out$gene_dist), c("Pr_0", "Pr_1", "Pr_2")),
                   out$par)
      snpdf <- rbind(snpdf, as.data.frame(snpdf_i))

      inddf_i <- data.frame(snp = snp_names[i], ind = ind_names, ref = out$input$refvec, size = out$input$sizevec,
                            geno = out$geno, postmean = out$postmean, maxpostprob = out$maxpostprob, row.names = NULL)
      inddf_i <- cbind(inddf_i, `colnames<-`(out$postmat, c("Pr_0", "Pr_1", "Pr_2")), `colnames<-`(out$genologlike, c("logL_0", "logL_1", "logL_2")))
      inddf <- rbind(inddf, inddf_i)
    }
    # Return
    return(list(snpdf = snpdf, inddf = inddf))
  }

  cat("\nRunning updog using", n.core, "cores...")
  if (n.core == 1) {
    out <- multidog(refmat = refmat, sizemat = sizemat, ploidy = 2, model = "norm")

  } else {
    # Split the markers into cores
    idx <- seq_along(mar_keep2)
    snp_cors <- split(x = idx, f = cut(x = idx, breaks = n.core))
    matlist <- lapply(X = snp_cors, FUN = function(ii) list(rmat = refmat[ii, ], smat = sizemat[ii, ]))
    # make clusters and apply a function over this list
    clust <- makeCluster(n.core)
    # Apply
    out <- parLapply(clust, matlist, run_flexdog)
    # out <- pblapply(X = matlist, FUN = run_flexdog, cl = clust)
    stopCluster(clust)

  }

  cat("done!")


  flexdog_list <- list()
  for (i in 1:100) {
    out <- flexdog(refvec = refmat1[i, ], sizevec = sizemat1[i, ], ploidy = 2, model = "norm")
  }
  multidog_out <- multidog(refmat = refmat[1:10,],sizemat = sizemat[1:10,], ploidy = 2, model = "norm", nc = 1)



}







