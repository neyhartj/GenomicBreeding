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
