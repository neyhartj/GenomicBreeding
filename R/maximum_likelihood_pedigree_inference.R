# mendel_prob <- function(g_child, g_mom, g_dad) {
#   # Convert genotype dosage (0,1,2) to allele counts
#   # Return probability of child genotype given parents
#
#   # All possible gametes from each parent
#   gametes_mom <- c(rep(0, 2 - g_mom), rep(1, g_mom))
#   gametes_dad <- c(rep(0, 2 - g_dad), rep(1, g_dad))
#
#   # All possible zygotes
#   zygotes <- outer(gametes_mom, gametes_dad, "+")
#
#   # Probability = fraction of zygotes matching child genotype
#   mean(zygotes == g_child)
# }
#
# # Log-likelihood for one individual
# loglik_individual <- function(geno_i, geno_mom, geno_dad, allele_freqs) {
#   ll <- 0
#   if (any(is.na(allele_freqs))) stop("Cannot have missing allele frequencies")
#
#   if (all(is.na(geno_mom)) || all(is.na(geno_dad))) {
#     # Founder: use Hardy-Weinberg
#     p <- allele_freqs
#
#     probs <- cbind("1" = (1-p)^2, "2" = 2*p*(1-p), "3" = p^2)
#     idx <- cbind(which(!is.na(geno_i)), geno_i[which(!is.na(geno_i))])
#
#     ll <- sum(log(probs[idx] + 1e-12))
#
#   } else {
#     for (l in which(!is.na(geno_i))) {
#       # Non-founder: Mendelian transmission
#       prob <- mendel_prob(g, geno_mom[l], geno_dad[l])
#       ll <- ll + log(prob + 1e-12)  # avoid log(0)
#
#     }
#
#   }
#
#   ll
# }
#
# # Log-likelihood for an entire pedigree
# loglik_pedigree <- function(ped, geno, allele_freqs) {
#   ll <- 0
#
#   for (i in seq_len(nrow(ped))) {
#     id  <- ped$id[i]
#     mom <- ped$mom[i]
#     dad <- ped$dad[i]
#
#     geno_i <- geno[id, ]
#     idx <- which(!is.na(geno_i))
#     if (!is.na(mom)) {
#       geno_mom <- geno[mom, ]
#       idx <- intersect(idx, which(!is.na(geno_mom)))
#     } else {
#       geno_mom <- rep(NA, length(geno_i))
#     }
#     if (!is.na(dad)) {
#       geno_dad <- geno[dad, ]
#       idx <- intersect(idx, which(!is.na(geno_dad)))
#     } else {
#       geno_dad <- rep(NA, length(geno_i))
#     }
#     geno_i <- geno_i[idx]
#     geno_dad <- geno_dad[idx]
#     geno_mom <- geno_mom[idx]
#
#     ll <- ll + loglik_individual(geno_i, geno_mom, geno_dad, allele_freqs)
#   }
#
#   ll
# }
#
# find_ml_pedigree <- function(pedigrees, geno, allele_freqs) {
#   ll <- sapply(pedigrees, loglik_pedigree, geno = geno, allele_freqs = allele_freqs)
#   best <- pedigrees[[which.max(ll)]]
#   list(best_pedigree = best, loglik = ll)
# }
#
#
# loglik_GA <- function(G, A, sigma2 = 0.01) {
#   diff <- G - A
#   ss <- sum(diff * diff)
#   -0.5 * ss / sigma2
# }
#
# set.seed(123)
#
# n_founders <- 50
# n_gen1     <- 25   # first generation after founders
# n_gen2     <- 25   # second generation
#
# # Founder IDs
# founders <- as.character(seq_len(n_founders))
#
# # Generation 1 IDs
# gen1 <- as.character(seq(n_founders + 1, n_founders + n_gen1))
#
# # Generation 2 IDs
# gen2 <- as.character(seq(n_founders + n_gen1 + 1, n_founders + n_gen1 + n_gen2))
#
# # Gen1 parents sampled from founders
# parents_gen1 <- t(replicate(n_gen1, sample(founders, 2, replace = FALSE)))
#
# # Gen2 parents sampled from Gen1
# parents_gen2 <- t(replicate(n_gen2, sample(gen1, 2, replace = FALSE)))
#
# # Build pedigree
# bigPed <- data.frame(
#   id   = c(founders, gen1, gen2),
#   dam  = c(rep(NA, n_founders), parents_gen1[,1], parents_gen2[,1]),
#   sire = c(rep(NA, n_founders), parents_gen1[,2], parents_gen2[,2]),
#   stringsAsFactors = FALSE
# )
#
# bigPed
#
# bigPed_df <- data.frame(sire = bigPed$sire, dam = bigPed$dam, label = bigPed$id)
#
# bigPed$sire[bigPed$id == "85"] <- bigPed$dam[bigPed$id == "85"]
#
# bigPed
#
# geno <- SimGeno(
#   Pedigree = bigPed,
#   nSnp = 200,
# )
# geno[geno == -9] <- as.numeric(NA)
#
# founders <- subset(bigPed, is.na(dam) & is.na(sire), id, drop = TRUE)
# allele_freqs <- colMeans(geno[founders,], na.rm = TRUE) / 2
# # G <- polyBreedR::G_mat(geno = t(geno), ploidy = 2, p.ref = allele_freqs)
# G <- polyBreedR::G_mat(geno = t(geno), ploidy = 2)
#
#
# ped_full <- pedigree(sire = bigPed$sire, dam = bigPed$dam, label = bigPed$id)
# A <- getA(ped_full)
#
# label <- "85"
# ped_cross <- bigPed_df[bigPed_df$label == label, ]
# ped_cross <- rbind(bigPed_df[bigPed_df$label %in% unlist(ped_cross[c(1, 2)]),], ped_cross)
# ped_cross$sire[1:2] <- as.character(NA)
# ped_cross$dam[1:2] <- as.character(NA)
# ped_cross
#
# ped_cross1 <- pedigree(sire = ped_cross$sire, dam = ped_cross$dam, label = ped_cross$label)
# Across <- getA(ped_cross1)
# ped_unre <- ped_cross
# ped_unre$sire <- as.character(NA)
# ped_unre$dam <- as.character(NA)
# ped_unre <- pedigree(sire = ped_unre$sire, dam = ped_unre$dam, label = ped_unre$label)
# Aunre <- getA(ped_unre)
# ped_self <- ped_cross
# ped_self$sire[3] <- ped_self$dam[3]
# ped_self <- pedigree(sire = ped_self$sire, dam = ped_self$dam, label = ped_self$label)
# ped_self
# Aself <- getA(ped_self)
#
# Alist <- list(cross = Across, unrelated = Aunre, self = Aself)
#
# sapply(Alist, FUN = function(Amat) loglik_GA(G = G[row.names(Amat), row.names(Amat)], A = Amat))
#
#
#
#
#
#
# ped_true <- data.frame(id = c("8", "13", "14"), mom = c(NA, NA, "8"), dad = c(NA, NA, "8"))
# ped_unr <- data.frame(id = c("8", "13", "14"), mom = c(NA, NA, NA), dad = c(NA, NA, NA))
# ped_alt1 <- data.frame(id = c("8", "13", "14"), mom = c(NA, NA, "8"), dad = c(NA, NA, "13"))
#
# ped_true <- data.frame(id = c("9", "10", "13"), mom = c(NA, NA, "9"), dad = c(NA, NA, "10"))
# ped_unr <- data.frame(id = c("9", "10", "13"), mom = c(NA, NA, NA), dad = c(NA, NA, NA))
# ped_alt1 <- data.frame(id = c("9", "10", "13"), mom = c(NA, NA, "9"), dad = c(NA, NA, "9"))
#
#
# ped_list <- list(true = ped_true, unrelated = ped_unr, alt1 = ped_alt1)
#
#
# find_ml_pedigree(pedigrees = ped_list, geno, allele_freqs)
