#' Compute (and plot) linkage disequilibrium decay
#'
#' @param geno A data frame of genotype calls, where the first three columns are "marker," "chrom,", and "position."
#' Subsequent columns are genotype calls for individuals.
#' @param r2.df A data frame of pairwise LD estimates, where the first four columns are "chrom", "marker1", "marker2",
#' and "r2". This can be obtained from the output of \code{calc_LD_r2}.
#' @param max.pair The maximum number of marker pairs to use when fitting the LD decay model.
#' @param xlim The limits to the x-axis when plotting LD decay.
#' @param method The method of computing pairwise LD (r2). If \code{"full"}, complete pairwise correlation coefficients
#' are computed; if \code{"approximate"}, an VERY approximate matrix multiplication method is used.
#'
#' @importFrom broom tidy
#' @importFrom dplyr bind_rows
#' @import ggplot2
#'
#' @export
#'
calc_LD_r2 <- function(geno, r2.df = NULL, by.chrom = TRUE, max.pair = 100000, xlim = c(0, 1e5),
                       method = c("full", "approximate")) {

  stopifnot(is.data.frame(geno))
  cols.match <- names(geno)[1:3] == c("marker", "chrom", "position")
  if (any(!cols.match)) {
    stop("The first 3 columns of 'geno' should be 'marker', 'chrom', and 'position'.")
  }

  if (!is.null(r2.df)) {
    stopifnot(is.data.frame(r2.df))
    cols.match <- names(r2.df)[1:4] == c("chrom", "marker1", "marker2", "r2")
    if (any(!cols.match)) {
      stop("The first four columns of 'r2.df' should be 'chrom', 'marker1', 'marker2', and 'r2'")
    }
  }

  stopifnot(is.logical(by.chrom))
  stopifnot(length(xlim) == 2)
  method <- match.arg(method)

  snp_info_use <- geno[,c("marker", "chrom", "position")]


  # If r2.df is provided, do not compute LD
  if (is.null(r2.df)) {

    # Convert to a matrix
    geno_mat <- rrblup2genomat(x = geno)

    # Remove markers with no variance
    geno_mat_sd <- apply(X = geno_mat, MARGIN = 2, FUN = sd, na.rm = TRUE)
    geno_mat <- geno_mat[, geno_mat_sd > 0, drop = FALSE]

    # Split markers by chromosome
    snp_info_use <- snp_info_use[snp_info_use$marker %in% colnames(geno_mat), ]
    markers_by_chrom <- split(snp_info_use, snp_info_use$chrom)

    # A list of tidy correlation matrices
    tidy_ld_r2 <- list()
    # Iterate over chromosomes
    for (j in seq_along(markers_by_chrom)) {
      markers_j <- markers_by_chrom[[j]]$marker
      geno_mat_j <- geno_mat[, markers_j]

      # Compute the correlations
      if (method == "full") {
        r2_ld_j <- cor(x = geno_mat_j, use = "pairwise.complete.obs")^2

      } else {
        # Fill in missing data with the mean
        geno_mat_j <- apply(X = geno_mat_j, MARGIN = 2, FUN = function(snp){
          snp[is.na(snp)] <- mean(snp, na.rm = TRUE)
          snp
        })

        geno_mat_j <- scale(geno_mat_j)
        r2_ld_j2 <- (crossprod(geno_mat_j) / (nrow(geno_mat_j) - 1))^2

      }

      # Tidy up the matrix
      r2_ld_j2 <- as.dist(r2_ld_j2)
      r2_ld_j2 <- broom::tidy(r2_ld_j2)
      names(r2_ld_j2) <- c("marker1", "marker2", "r2")
      r2_ld_j2 <- cbind(chrom = names(markers_by_chrom)[[j]], r2_ld_j2)

      tidy_ld_r2[[names(markers_by_chrom)[[j]]]] <- r2_ld_j2

    }

    # Merge and return
    r2.df <- bind_rows(tidy_ld_r2)

  }

  # Distance between markers
  snp_info_use <- snp_info_use[snp_info_use$marker %in% c(r2.df$marker1, r2.df$marker2), ]
  snp_pos <- as.matrix(snp_info_use$position)
  row.names(snp_pos) <- snp_info_use$marker
  r2.df$bp_dist <- abs(snp_pos[r2.df$marker1, ] - snp_pos[r2.df$marker2, ])


  ## Fit the non-linear models of Abecasis et al. 2001 as coded by Anderson et al 2016:
  start <- list(dlow=0.05, dhigh=0.95, gens=5000)
  lower <- c(0.001, 0.001, 1)
  upper <- c(0.999, 0.999, 100000)

  # Use a recombination rate of 1 cM / 1 Mb
  recomb_rate <- 1
  # Scaling factor
  scaling <- (1e6 * 100 * 10)
  # Distance for plotting
  dist_max <- max(xlim)

  # Split by chromosome or no?
  if (by.chrom) {
    r2_df_split <- split(r2.df, r2.df$chrom)

  } else {
    r2_df_split <- list(r2.df)

  }

  ## Model LD to create the curves
  r2_model_list <- list()
  # Iterate over the split
  for (j in seq_along(r2_df_split)) {
    df <- r2_df_split[[j]]

    # Randomly sample pairs
    max.pair.j <- ifelse(max.pair > nrow(df), nrow(df), max.pair)
    df <- df[sort(sample(seq_len(nrow(df)), max.pair.j)), ]

    ld <- df$r2
    d <- df$bp_dist

    # Calculate recombination fraction (Mb * (cM/Mb)) / 100 cM;
    # Use this same formula in the model for r
    r <- (((d / 1e6) * recomb_rate) / 100)

    # Fit the model
    ld_decay_model <- nls(formula = ld ~ dlow + (dhigh - dlow) * ( (1 - (d / scaling) * 1)^gens),
                          start = start, control = nls.control(maxiter = 1000, warnOnly = T),
                          lower = lower, upper = upper, algorithm = "port", trace = FALSE)

    # Predict
    xvals <- seq(0, 1e7, 500)
    yhat <- predict(object = ld_decay_model, list(d = xvals))
    predictions <- data.frame(d = xvals, d_kbp = xvals / 1000,  yhat = yhat)

    # ## plot?
    # plot(ld ~ d)
    # points()

    model_summary <- summary(ld_decay_model)
    dlow_out <- model_summary$parameters["dlow", "Estimate"]
    dhigh_out <- model_summary$parameters["dhigh", "Estimate"]
    gens_out <- model_summary$parameters["gens", "Estimate"]

    # Calculate the half-life
    half <- yhat[1]/2
    half_scaled <- 1 - 10^(log10((half - dlow_out)/(dhigh_out - dlow_out))/gens_out)
    half_bp <- (half_scaled * scaling) / recomb_rate

    r2_model_list[[j]] <- list(predictions = cbind(chrom = names(r2_df_split)[j], predictions), half.life.bp = half_bp)

  }


  # Plot!
  r2_model_df <- do.call(rbind, lapply(r2_model_list, "[[", "predictions"))
  half_bp <- sapply(r2_model_list, "[[", "half.life.bp")

  g <- ggplot(data = r2_model_df, aes(x = d_kbp, y = yhat))
  if (by.chrom) {
    g <- g + geom_line(aes(color = chrom))
  } else {
    g <- g + geom_line()
  }
  g <- g +
    scale_x_continuous(name = "Distance (kbp)", breaks = pretty, limits = xlim / 1000) +
    scale_y_continuous(name = expression("Linkage disequilibrium ("*italic(r^2)*")"), breaks = pretty) +
    ggsci::scale_color_igv(name = "Chromosome") +
    theme_classic() +
    theme(legend.position = "inside", legend.position.inside = c(0.8, 0.9), legend.justification = c(0, 1))


  cat("\n\nSummary of LD decay half-life (bp):\n")
  print(summary(half_bp))

  ld_decay_thresh <- c(0.3, 0.2, 0.1)
  names(ld_decay_thresh) <- ld_decay_thresh
  ld_decay_thresh <- sapply(X = ld_decay_thresh, FUN = function(threshold) {
    r2_below_thresh <- r2_model_df[r2_model_df$yhat <= threshold, ]
    tapply(X = r2_below_thresh$d, INDEX = r2_below_thresh$chrom, min)
  })

  ld_decay_thresh_summ <- colMeans(ld_decay_thresh)
  ld_decay_thresh_summ <- data.frame(LD.threshold = colnames(ld_decay_thresh), bp = round(ld_decay_thresh_summ), row.names = NULL)

  cat("\n\nSummary of average bp distance to reach LD levels:\n")
  print(ld_decay_thresh_summ)

  # Plot the chart
  print(g)

  half_bp_out <- matrix(half_bp, dimnames = list(row.names(ld_decay_thresh), "bp"))

  # Return
  output <- list(r2.df = r2.df, decay.predictions = r2_model_df, half.life.bp = half_bp_out,
                 ld.decay.thresholds = ld_decay_thresh, plot = g)

  return(output)

}

