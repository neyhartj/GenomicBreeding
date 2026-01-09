#' Prep a qtl2 cross object for QTL modelling
#'
#'
#'
#'
prep_cross2 <- function(cross, ped, ploidy = 2, error.prob = 0.05, cores = 1, crosstype = c("f2", "outcross"),
                        probs = NULL, map = NULL) {

  stopifnot(inherits(cross, "cross2"))
  stopifnot(ploidy %in% c(2, 4))
  stopifnot(is.numeric(error.prob))
  stopifnot(cores >= 1)
  stopifnot(is.data.frame(ped))
  crosstype <- match.arg(crosstype)

  # If probs is NULL, calculate genotype probabilities
  if (is.null(probs)) {

    # Calculate genotype probabilities
    cat("\nCalculating genotype probabilities...")
    probs <- calc_genoprob(cross = cross, map = cross$gmap, error_prob = error.prob, cores = cores)
    alleleprobs <- genoprob_to_alleleprob(probs = probs)


  } else {
    # Convert genotype probs to allele probs if needed
    stopifnot(inherits(probs, "calc_genoprob"))
    if (!attr(probs, "alleleprobs")) {
      alleleprobs <- genoprob_to_alleleprob(probs = probs)
    } else {
      alleleprobs <- probs
    }

  }

  # Get map information from the cross
  if (is.null(map)) {
    # get the marker map information
    map_info <- qtl::map2table(cross$gmap)
    map_info <- cbind(marker = row.names(map_info), map_info)
    names(map_info) <- c("marker", "chrom", "cM")
    # Add physical map information if available
    if (!is.null(cross$pmap)) {
      pmap_info <- qtl::map2table(map = cross$pmap)
      pmap_info <- cbind(marker = row.names(pmap_info), pmap_info)
      map_info$bp <- pmap_info$pos

    }

    map <- map_info

  } else {
    stopifnot(names(map)[1:3] == c("marker", "chrom", "cM"))

  }

  # Total number of markers
  nmar <- sum(n_mar(cross))

  # Compute the additive dosage matrix. This is alleleprob * ploidy
  geno_add <- lapply(X = alleleprobs, FUN = function(x) {
    chr_list <- list()
    marnames <- dimnames(x)[[3]]
    for (k in seq_along(marnames)) {
      chr_list[[marnames[k]]] <- x[,,k] * ploidy
    }
    chr_list
  })
  geno_add <- do.call(c, geno_add)
  names(geno_add) <- marker_names(cross)

  # Compute the dominance dosage matrix.
  #
  # if crosstype is F2 and only allele probs were given, stop
  if (crosstype == "f2" & attr(probs, "alleleprobs")) {
    stop("You must provide genotype probabilities for an F2.")
  }

  # For f2s, use the genotype probabilites.
  if (crosstype == "f2") {
    geno_dom <- lapply(X = probs, FUN = function(x) {
      chr_list <- list()
      marnames <- dimnames(x)[[3]]
      for (k in seq_along(marnames)) {
        chr_list[[marnames[k]]] <- as.matrix(x[,2, k])
      }
      chr_list
    })

    genoprobs <- probs

  } else {
    stop("Support for outcrosses is not available yet.")
  }

  geno_dom <- do.call(c, geno_dom)
  names(geno_dom) <- marker_names(cross)

  # Create the geno object
  geno <- vector("list", length = length(geno_dom))
  for (k in seq_along(geno)) {
    geno[[k]] <- list(geno_add[[k]], geno_dom[[k]])
  }
  names(geno) <- names(geno_add)

  # Create a list to output
  output <- list(
    X = matrix(1, nrow = nrow(probs[[1]]), ncol = 1, dimnames = list(NULL, "(Intercept)")),
    Z = diag(nrow(probs[[1]])),
    map = map,
    geno = geno,
    probs = genoprobs,
    crosstype = crosstype
  )

  return(output)

}



# # Test objects
# cross <- cross_qtl2_STxOxy
#
# trait <- "DaysToFirstFlw"
# probs <- calc_genoprob(cross = cross, map = cross$gmap, error_prob = 0.06, cores = 12)
# kin <- calc_kinship(probs = probs, type = "loco", cores = 12)
# scan_out <- qtl2::scan1(genoprobs = probs, pheno = cross$pheno[,trait, drop = FALSE], kinship = kin)
# peaks <- find_peaks(scan1_output = scan_out, map = cross$pmap, threshold = 4)
# qtl_mars <- find_marker(map = cross$pmap, chr = peaks$chr, pos = peaks$pos)
#
# ploidy <- 2
# qtl <- data.frame(marker = qtl_mars, dominance = 2)
# epistasis <- as.data.frame(t(combn(x = qtl_mars, 2)))
# names(epistasis) <- c("marker1", "marker2")
# polygenic <- TRUE
# ped <- data.frame(id = ind_ids(cross), parent1 = "STEVENS", parent2 = "NJ96-20")
#
# # Convert the cross object
# cross_diaqtl <- prep_cross2(cross = cross, ped = ped, ploidy = 2, error.prob = 0.06, cores = 12)
#
#
#
#
# Fit the QTL model
fit_qtl_lmm <- function(data, pheno, trait, ploidy = 2, qtl, epistasis = NULL, polygenic = TRUE, max.cor.sing = 0.97) {


  # error checking
  stopifnot(is.character(trait))
  stopifnot(is.data.frame(qtl))
  stopifnot(is.logical(polygenic))
  stopifnot(ploidy %in% c(2, 4))
  stopifnot(is.matrix(pheno))

  # Order of genos and phenos is assumed to be the same


  # Function for formatting epistasis
  format_Xaa <- function(j, k, probs.list) {
    states <- colnames(probs[[1]])
    states <- as.integer(as.factor(states))
    genoprob1 <- probs.list[[j]]
    genoprob2 <- probs.list[[k]]

    # This is hard-coded for F2s
    u1 <- cbind(c(2, 0), c(1, 1), c(0, 2))
    uu <- kronecker(u1, u1)

    tmp3 <- sapply(seq_len(nrow(genoprob1)), FUN = function(i) {
      w1 <- genoprob1[i, ]
      w2 <- genoprob2[i, ]
      ww <- matrix(kronecker(w1, w2), ncol = 1)
      uu %*% ww
    })

    return(t(tmp3))

  }

  # Make sure the markers are listed
  all_markers <- data$map$marker
  stopifnot(qtl$marker %in% all_markers)
  # Dominance cannot exceed ploidy
  if (max(qtl$dominance) > ploidy) stop("Dominance levels cannot be higher than the ploidy.")

  markers <- unique(qtl$marker)
  n.qtl <- length(markers)

  y <- pheno[, trait]
  if (inherits(y, "factor")) {
    response <- "ordinal"
  } else {
    response <- "gaussian"
  }

  idx <- which(!is.na(y))
  y1 <- y[idx]

  Xqtl <- vector("list", n.qtl)
  for (i in seq_len(n.qtl)) {
    dominance <- qtl$dominance[i]
    Xqtl[[i]] <- vector("list", dominance)
    for (j in 1:dominance) {
      Xqtl[[i]][[j]] <- data$geno[[markers[i]]][[j]]
    }
  }

  # Format epistasis input
  if (!is.null(epistasis)) {
    # bind all genoprobs together
    probs_all <- lapply(X = probs, FUN = function(x) {
      chr_list <- vector("list", dim(x)[3])
      for (k in seq_along(chr_list)) {
        chr_list[[k]] <- x[,,k]
      }
      names(chr_list) <- dimnames(x)[[3]]
      chr_list
    })
    marnames_list <- lapply(probs_all, names)
    probs_all <- do.call(c, probs_all)
    names(probs_all) <- unlist(marnames_list, use.names = FALSE)

    marker1 <- epistasis$marker1
    marker2 <- epistasis$marker2
    if (length(setdiff(c(marker1, marker2), markers)) > 0) {
      stop("Markers in epistasis should also be in qtl")
    }
    n.epi <- nrow(epistasis)
    Xaa <- vector("list", n.epi)
    for (i in 1:n.epi) {
      Xaa[[i]] <- format_Xaa(j = marker1[i], k = marker2[i], probs.list = probs_all)
    }
  } else {
    n.epi <- 0
    marker1 <- marker2 <- Xaa <- NULL
  }

  if (polygenic) {
    poly.marker <- unique(c(markers, marker1, marker2))
    chroms <- setdiff(unique(data$map$chrom), data$map$chrom[match(poly.marker, data$map$marker)])
    if (length(chroms) == 0) {
      stop("There are no chromosomes remaining for the polygenic effect.")
    }
    polyG <- calc_kinship(probs = qtl2:::subset.calc_genoprob(x = probs, chr = chroms), type = "overall")
  } else {
    polyG <- diag(nrow(probs[[1]]))
  }

  polyG1 <- polyG[idx, idx]

  Xqtl1 <- Xqtl
  for (i in seq_along(Xqtl)) {
    x <- Xqtl[[i]]
    xa <- x[[1]][idx, , drop = FALSE]
    colnames(xa) <- paste0(markers[i], ".", colnames(xa))
    xd <-  x[[2]][idx, , drop = FALSE]
    colnames(xd) <- paste0(markers[i], ".dom")
    Xqtl1[[i]] <- list(xa, xd)
  }

  Xaa1 <- Xaa
  for (i in seq_along(Xaa)) {
    aa_i <- Xaa[[i]][idx, , drop = FALSE]
    colnames(aa_i) <- c(paste0(epistasis$marker1[i], ".1:", epistasis$marker2[i], ".1"),
                        paste0(epistasis$marker1[i], ".1:", epistasis$marker2[i], ".2"),
                        paste0(epistasis$marker1[i], ".2:", epistasis$marker2[i], ".1"),
                        paste0(epistasis$marker1[i], ".2:", epistasis$marker2[i], ".2"))
    Xaa1[[i]] <- aa_i
  }

  # Create a data frame for fitting the model
  dat <- data.frame(y = y1, id = row.names(polyG1))

  # Random formula for additive effects
  add_terms_all <- add_terms <- paste0("Xqtl1[[", seq_along(Xqtl1), "]][[1]]")
  dom_terms_all <- dom_terms <- paste0("Xqtl1[[", seq_along(Xqtl1), "]][[2]]")
  int_terms_all <- int_terms <- paste0("Xaa1[[", seq_along(Xaa1), "]]")
  rand <- reformulate(termlabels = c(add_terms, dom_terms, int_terms))

  poly_term <- "vsm(ism(id), Gu = polyG1)"

  # Calculate the corelation between additive model matrices
  add_mat_cor <- matrix(0, nrow = length(Xqtl1), ncol = length(Xqtl1), dimnames = list(qtl$marker, qtl$marker))
  for (i in seq_len(nrow(add_mat_cor))) {
    for (j in seq_len(ncol(add_mat_cor))) {
      add_mat_cor[i,j] <- mean(abs(cor(Xqtl1[[i]][[1]], Xqtl1[[j]][[1]])))
    }
  }
  # Remove pairs with > max.cor.sing correlation
  add_mat_cor_rm <- which(add_mat_cor > max.cor.sing, arr=TRUE)
  add_mat_cor_rm <- add_mat_cor_rm[add_mat_cor_rm[,1] < add_mat_cor_rm[,2],]
  if (nrow(add_mat_cor_rm) > 0) {
    add_terms_rm <- unique(add_mat_cor_rm[,2])
    add_terms <- add_terms[-add_terms_rm]
    # Print a message
    warning(paste0(length(add_terms_rm), " additive terms were removed to avoid singularities."))
  }

  # Repeat with dominance effects
  # Calculate the corelation between dominance model matrices
  dom_mat_cor <- matrix(0, nrow = length(Xqtl1), ncol = length(Xqtl1), dimnames = list(qtl$marker, qtl$marker))
  for (i in seq_len(nrow(dom_mat_cor))) {
    for (j in seq_len(ncol(dom_mat_cor))) {
      dom_mat_cor[i,j] <- mean(abs(cor(Xqtl1[[i]][[2]], Xqtl1[[j]][[2]])))
    }
  }
  # Remove pairs with > max.cor.sing correlation
  dom_mat_cor_rm <- which(dom_mat_cor > max.cor.sing, arr=TRUE)
  dom_mat_cor_rm <- dom_mat_cor_rm[dom_mat_cor_rm[,1] < dom_mat_cor_rm[,2],]

  # correlate dom terms with additive
  # Calculate the corelation between dominance model matrices
  add_dom_mat_cor <- matrix(0, nrow = length(Xqtl1), ncol = length(Xqtl1), dimnames = list(qtl$marker, qtl$marker))
  for (i in seq_len(nrow(add_dom_mat_cor))) {
    for (j in seq_len(ncol(add_dom_mat_cor))) {
      add_dom_mat_cor[i,j] <- mean(abs(cor(Xqtl1[[i]][[1]], Xqtl1[[j]][[2]])))
    }
  }
  # Remove pairs with > max.cor.sing correlation
  add_dom_mat_cor_rm <- which(add_dom_mat_cor > max.cor.sing, arr=TRUE)
  if (nrow(dom_mat_cor_rm) > 0 | nrow(add_mat_cor_rm) > 0 | nrow(add_dom_mat_cor_rm) > 0) {
    dom_terms_rm <- unique(c(dom_mat_cor_rm[,2], add_mat_cor_rm[,2], add_dom_mat_cor_rm[,2])) # remove dominance terms for additive terms that were removed.
    dom_terms <- dom_terms[-dom_terms_rm]
    # Print a message
    warning(paste0(length(dom_terms_rm), " dominance terms were removed to avoid singularities."))
  }

  # Repeat with epistatic effects
  # Calculate the corelation between epistatic model matrices
  epi_mat_cor <- matrix(0, nrow = length(Xaa1), ncol = length(Xaa1), dimnames = list(apply(epistasis, 1, paste0, collapse = "+"), apply(epistasis, 1, paste0, collapse = "+")))
  for (i in seq_len(nrow(epi_mat_cor))) {
    for (j in seq_len(ncol(epi_mat_cor))) {
      epi_mat_cor[i,j] <- mean(abs(diag(cor(Xaa1[[i]], Xaa1[[j]]))))
      # epi_mat_cor[i,j] <- mean(abs(cor(Xaa1[[i]], Xaa1[[j]])))

    }
  }
  # Remove pairs with > max.cor.sing correlation
  epi_mat_cor_rm <- which(epi_mat_cor > max.cor.sing, arr=TRUE)
  epi_mat_cor_rm <- epi_mat_cor_rm[epi_mat_cor_rm[,1] < epi_mat_cor_rm[,2],]
  if (nrow(epi_mat_cor_rm) > 0 | nrow(add_mat_cor_rm) > 0) {
    # Find the epistatic terms corresponding to the additive terms that were removed
    add_mar_rm <- names(add_mat_cor_rm[,2])
    epi_terms_add_mar_rm <- which(apply(sapply(epistasis, function(x) x %in% add_mar_rm), 1, any))
    # Now find epistatic terms to be removed
    epi_terms_rm <- unique(epi_mat_cor_rm[,2])
    # Combine
    epi_terms_rm <- union(epi_terms_add_mar_rm, epi_terms_rm)
    int_terms <- int_terms[-epi_terms_rm]
    # Print a message
    warning(paste0(length(epi_terms_rm), " epistatic terms were removed to avoid singularities."))
  }

  try_fit <- try(ans0 <- mmes(fixed = y ~ 1, random = reformulate(c(add_terms, dom_terms, int_terms, poly_term)), data = dat, verbose = FALSE), silent = TRUE)
  if (class(try_fit) == "try-error") stop("Singularities in the mixed-model. Try reducing 'max.cor.sing'.")

  # Start backwards elimination
  #
  # Epistatic effects first
  int_terms_test <- int_terms
  reduce <- TRUE
  while (reduce) {
    epi_terms_test_p <- numeric(length = length(int_terms_test))
    for (i in seq_along(int_terms_test)) {
      int_terms_dropped <- int_terms_test[-i]
      try_ans <- try(ans1 <- mmes(fixed = y ~ 1, random = reformulate(c(add_terms, dom_terms, int_terms_dropped, poly_term)), data = dat, verbose = FALSE), silent = TRUE)
      if (class(try_ans) == "try-error") {
        epi_terms_test_p[i] <- as.numeric(NA)
      } else {
        out <- capture.output(test <- anova(ans1, ans0))
        p <- as.numeric(sub(pattern = "\\*.*", replacement = "", x = test$PrChisq[2]))
        epi_terms_test_p[i] <- p
      }
    }
    if (any(epi_terms_test_p > 0.05)) {
      reduce <- TRUE
      int_terms_test <- int_terms_test[-which.max(epi_terms_test_p)]
    } else {
      reduce <- FALSE
      int_terms_use <- int_terms_test
    }
  }

  # Refit the model
  ans0.1 <- mmes(fixed = y ~ 1, random = reformulate(c(add_terms, dom_terms, int_terms_use, poly_term)), data = dat, verbose = FALSE)

  # dominance effects next
  dom_terms_test <- dom_terms
  reduce <- TRUE
  while (reduce) {
    terms_test_p <- numeric(length = length(dom_terms_test))
    for (i in seq_along(dom_terms_test)) {
      dom_terms_dropped <- dom_terms_test[-i]
      try_ans <- try(ans1 <- mmes(fixed = y ~ 1, random = reformulate(c(add_terms, dom_terms_dropped, int_terms_use, poly_term)), data = dat, verbose = FALSE), silent = TRUE)
      if (class(try_ans) == "try-error") {
        terms_test_p[i] <- as.numeric(NA)
      } else {
        out <- capture.output(test <- anova(ans1, ans0.1))
        out <- sub(pattern = "\\*.*$", replacement = "", x = test$PrChisq[2])
        out <- sub(pattern = "\\.$", replacement = "", x = out)
        p <- as.numeric(out)
        terms_test_p[i] <- p
      }
    }
    if (any(terms_test_p > 0.05)) {
      reduce <- TRUE
      dom_terms_test <- dom_terms_test[-which.max(terms_test_p)]
    } else {
      reduce <- FALSE
      dom_terms_use <- dom_terms_test
    }
  }

  ans0.2 <- mmes(fixed = y ~ 1, random = reformulate(c(add_terms, dom_terms_use, int_terms_use, poly_term)), data = dat, verbose = FALSE)

  # additive effects next
  # Do not test additive terms that are part of interactions or dominance
  which_dom_terms_use <- which(dom_terms_all %in% dom_terms_use)
  which_int_terms_use <- which(int_terms_all %in% int_terms_use)
  which_int_terms_use <- which(qtl$marker %in% unique(unlist(epistasis[which_int_terms_use, ])))
  which_add_terms_use <- which(add_terms_all %in% add_terms)

  add_terms_test <- add_terms_all[setdiff(which_add_terms_use, union(which_dom_terms_use, which_int_terms_use))]
  add_terms_keep <- add_terms_all[union(which_dom_terms_use, which_int_terms_use)]

  reduce <- TRUE
  while (reduce) {
    terms_test_p <- numeric(length = length(add_terms_test))
    for (i in seq_along(add_terms_test)) {
      add_terms_dropped <- sort(c(add_terms_keep, add_terms_test[-i]))
      try_ans <- try(ans1 <- mmes(fixed = y ~ 1, random = reformulate(c(add_terms_dropped, dom_terms_use, int_terms_use, poly_term)), data = dat, verbose = FALSE), silent = TRUE)
      if (class(try_ans) == "try-error") {
        terms_test_p[i] <- as.numeric(NA)
      } else {
        out <- capture.output(test <- anova(ans1, ans0.1))
        out <- sub(pattern = "\\*.*$", replacement = "", x = test$PrChisq[2])
        out <- sub(pattern = "\\.$", replacement = "", x = out)
        p <- as.numeric(out)
        terms_test_p[i] <- p
      }
    }
    if (any(terms_test_p > 0.05)) {
      reduce <- TRUE
      add_terms_test <- add_terms_test[-which.max(terms_test_p)]
    } else {
      reduce <- FALSE
      add_terms_use <- sort(c(add_terms_keep, add_terms_test))
    }
  }

  # Final model
  final_model <- mmes(fixed = y ~ 1, random = reformulate(c(add_terms_use, dom_terms_use, int_terms_use, poly_term)), data = dat, verbose = FALSE)
 # get effects
  u_hat <- final_model$u

  # Instantiate objects for the effects
  variances <- NULL
  haplotypes <- colnames(Xqtl[[1]][[1]])
  n.hap <- length(haplotypes)
  additive.effects <- matrix(as.numeric(NA), nrow = n.qtl, ncol = n.hap)
  colnames(additive.effects) <- haplotypes
  rownames(additive.effects) <- qtl$marker

  if (any(qtl$dominance > 1)) {
    diplotypes <- paste0(haplotypes, collapse = ".")
    n.diplo <- length(diplotypes)
    digenic.effects <- matrix(as.numeric(NA), nrow = n.qtl, ncol = n.diplo)
    colnames(digenic.effects) <- diplotypes
    rownames(digenic.effects) <- qtl$marker
  }

  for (i in 1:n.qtl) {
    qtl_list <- list()
    dominance <- qtl$dominance[i]
    tmp <- matrix(0, nrow = 1, ncol = dominance)
    colnames(tmp) <- paste(qtl$marker[i], c("additive", "digenic", "trigenic", "quadrigenic")[1:dominance])
    for (j in 1:dominance) {
      pattern <- switch(j,
        `1` = paste0(qtl$marker[i], ".", haplotypes),
        `2` = paste0(qtl$marker[i], ".dom")
      )
      idx_u <- match(x = pattern, table = row.names(u_hat))
      if (length(idx_u) == 0) {
        u_hat_j <- matrix(0, nrow = 1, ncol = 1)
      } else {
        u_hat_j <- t(u_hat[idx_u, ])
      }
      qtl_list[[j]] <- u_hat_j
      tmp[, j] <- var(tcrossprod(Xqtl[[i]][[j]], u_hat_j))
    }
    variances <- cbind(variances, tmp)

    additive.effects[i, ] <- qtl_list[[1]]
    if (dominance > 1) {
      digenic.effects[i, ] <- qtl_list[[2]]
    }
  }

  if (!is.null(epistasis)) {
    tmp <- matrix(0, nrow = 1, ncol = n.epi)
    colnames(tmp) <- paste(paste(epistasis$marker1, epistasis$marker2, sep = "+"), "epistasis")
    epistasis.effects <- matrix(NA, nrow = ncol(Xaa[[1]]), ncol = n.epi)
    for (i in 1:n.epi) {
      u_hat_epi <- t(u_hat[grep(pattern = paste0(epistasis$marker1[i], "\\..*:", epistasis$marker2[i]), row.names(u_hat)), ])
      if (ncol(u_hat_epi) == 0) {
        next
      } else {
        tmp[, i] <- var(tcrossprod(Xaa[[i]], u_hat_epi))
        epistasis.effects[, i] <- u_hat_epi
      }
    }
    colnames(epistasis.effects) = paste(paste(epistasis$marker1, epistasis$marker2, sep = "+"))
    rownames(epistasis.effects) = apply(expand.grid(haplotypes, haplotypes)[, 2:1], 1, paste0, collapse = "+")
    variances <- cbind(variances, tmp)
  }
  if (polygenic) {
    theta <- final_model$theta
    theta <- theta[grep(pattern = "polyG1", x = names(theta))]
    variances <- cbind(variances, polygenic = theta[[1]][1,1] * mean(diag(polyG)))
  }
  # Remove 0 or NA
  variances <- variances[,!(is.na(variances) | variances == 0), drop = FALSE]
  if (params$response == "gaussian") {
    variances <- cbind(variances, residual = tail(final_model$theta, 1)[[1]][1,1])
    h2 <- variances/rowSums(variances)
  } else {
    y2 <- as.integer(y) - 1
    R2 <- 1 - sum(ans1$resid^2, na.rm = T)/sum((y2 - mean(y2, na.rm = T))^2, na.rm = T)
    h2 <- variances/apply(variances, 1, sum) * R2
  }
  h2 <- t(h2)
  return.var <- h2

  if (max(qtl$dominance) > 1) {
    effects <- list(additive = t(additive.effects), digenic = t(digenic.effects))
  } else {
    effects <- list(additive = t(additive.effects))
  }
  if (!is.null(epistasis)) {
    effects[["epistasis"]] = epistasis.effects
  }
  return(list(resid = tail(final_model$theta, 1)[[1]][1,1], var = return.var,
              effects = effects))
}




#' Plot the results of fitted QTL model
#'
#' @param fit1out The output from \code{\link{fit1}[qtl]}.
#' @param genoprobs Genotype probabilities.
#'
#' @import ggplot2
#'
#' @export
#'
plot_fit1 <- function(fit1out, genoprobs) {
  coef_df <- data.frame(haplotype = names(fit1out$coef), effect = fit1out$coef, se = fit1out$SE, row.names = NULL)
  haps <- if (is.list(genoprobs)) {
    dimnames(genoprobs[[1]])[[2]]
  } else {
    dimnames(genoprobs)[[2]]
  }

  coef_df <- coef_df[coef_df$haplotype %in% haps, ]
  hap_split <- strsplit(x = coef_df$haplotype, split = "_hap", fixed = TRUE)
  coef_df$geno_id <- sapply(hap_split, "[[", 1)
  coef_df$haplotype <- sapply(hap_split, "[[", 2)
  coef_df$haplotype <- as.factor(coef_df$haplotype)

  g <- ggplot(data = coef_df, aes(x = haplotype, y = effect)) +
    geom_col(fill = "firebrick", color = "black") +
    geom_errorbar(aes(ymin = effect - se, ymax = effect + se), width = 0.3) +
    facet_grid(~ geno_id) +
    theme_bw()

  return(g)

}











