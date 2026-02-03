#' Principal component analysis of marker genotypes
#'
#'
#' @param geno A data.frame of genotype calls, where the first three columns are "marker," "chrom,", and "position."
#' Subsequent columns are genotype calls for individuals. The \code{geno} object is assumed to be filtered.
#' @param metadata Metadata for the individuals in \code{geno}. Used for plotting. The first column
#' of \code{metadata} must be "geno_id" for the individual name. Other columns are metadata columns.
#' @param input What should be the input for PCA? Can be \code{relmat} and an additive relationship
#' matrix will be computed, or \code{geno.mat} and PCA will be performed on the marker genotype matrix.
#' @param ploidy Integer ploidy number.
#' @param n.core The number of cores for processing.
#'
#'
#'
geno_PCA <- function(geno, metadata, input = c("relmat", "geno.mat"), ploidy = 2, n.core = 1) {

  stopifnot(is.data.frame(geno))
  cols.match <- names(geno)[1:3] == c("marker", "chrom", "position")
  if (any(!cols.match)) {
    stop("The first 3 columns of 'geno' should be 'marker', 'chrom', and 'position'.")
  }

  stopifnot(is.data.frame(metadata))
  if (names(metadata)[1] != "geno_id") stop("The first column of 'metadata' must be 'geno_id'.")

  stopifnot(ploidy %in% c(2, 4, 6))

  input <- match.arg(input)

  stopifnot(n.core < 1)
  n.core <- as.integer(n.core)

  # Convert the geno object to a matrix
  geno1 <- impute_geno(geno = geno, ploidy = ploidy, method = "pop", ind.max.missing = 1, snp.max.missing = 1, n.core = n.core)
  geno1 <- rrblup2genomat(x = geno1, transpose = FALSE)

  # Compute relationship matrix?
  if (input == "relmat") {
    pca_input <- G_mat(geno = geno1, ploidy = ploidy)
  } else {
    pca_input <- geno1
  }

  pca_ans <- prcomp(x = pca_input)

  pca_tidy <- pca_ans$x

  pca_tidy <- broom::tidy(x = K1_pca) %>%
    mutate(PC = paste0("PC", PC),
           PC = fct_inorder(PC)) %>%
    spread(PC, value) %>%
    rename(sample_id = row) %>%
    mutate(geno_id = sub(pattern = "\\.[0-9]{1,}$", replacement = "", x = sample_id)) %>%
    left_join(., germplasm) %>%
    mutate(origin_state = ifelse(is.na(origin_state), "Unknown", origin_state),
           origin_state = as.factor(origin_state),
           origin_state = fct_relevel(origin_state, "Unknown", after = Inf))


  # prop_var <- K1_pca$sdev^2
  # prop_var <- prop_var / sum(prop_var)
  # prop_var <- round(prop_var * 100, 3)
  # names(prop_var) <- paste0("PC", seq_along(prop_var))
  #
  # x <- "PC1"
  # y <- "PC2"
  #
  # g <- K1_pca_tidy %>%
  #   subset(category == "WILD") %>%
  #   ggplot(aes(x = .data[[x]], y = .data[[y]], color = origin_state, shape = organism, label = geno_id)) +
  #   geom_point() +
  #   scale_shape_discrete(guide = "none") +
  #   ggsci::scale_color_igv() +
  #   scale_x_continuous(name = paste0(x, " (", prop_var[x], "%)"), breaks = pretty, guide = guide_axis(cap = TRUE)) +
  #   scale_y_continuous(name = paste0(y, " (", prop_var[y], "%)"), breaks = pretty, guide = guide_axis(cap = TRUE)) +
  #   theme_classic()
  # g
  #
  # ggsave(filename = "DCran_DArT_wild_PCA_1.png", plot = g, path = getwd(), width = 8, height = 5, dpi = 300)
  #
  # plotly::ggplotly(p = g)

}
