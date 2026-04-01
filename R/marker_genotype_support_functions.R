#' Interpolate marker genetic positions
#' 
#' @description 
#' Uses a reference map of recombination rates (cM/Mb) to interpolate the genetic position of markers based on their physical positions.
#' 
#' @param data A data.frame of SNPs with positions to be extrapolated. This should resemble the table format of genetic maps from the R/qtl package: first column is chromosome, second column is physical position, row.names are marker names.
#' @param reference A data.frame of physical map intervals with smoothed recombination rates. Columns should be: 1) chromosome, 2) left physical position boundary, 3) right physical position boundary, and 4) recombination rate (in cM/Mb).
#' @param max.recomb.rate The maximum value of smoothed recombination rates. Any values in \code{reference} that exceed this value will be converted to this value.
#' 
#' @export
#' 
interpolate_positions <- function(x, reference, max.recomb.rate = 7) {
  
  ## Check formats
  stopifnot(is.data.frame(x))
  stopifnot(is.data.frame(reference))
  
  ## Replace anything greater than the max rate with the max
  reference[[4]][reference[[4]] > max.recomb.rate] <- max.recomb.rate
  
  ## Make sure the first columns of each input match
  if (!isTRUE(all.equal(unique(x[[1]]), unique(reference[[1]])))) 
    stop("The chromosome levels in 'x' and 'reference' do not match.")
  
  
  ## Split x and reference by chromosome
  x_split <- split(x, x[[1]])
  reference_split <- split(reference, reference[[1]])
  
  ## Simultaneous mapping
  interp_out <- mapply(x_split, reference_split, FUN = function(.x, .y) {
    
    # Add a 0 row
    .y1 <- rbind(.y[1,], .y) 
    .y1[1,2:3] <- c(0, .y1[[2]][1])
    
    # Add cm bounds to .y
    .y1$left_cM <- 0
    .y1$right_cM <- 0
    
    for (i in seq(nrow(.y1))) {
      
      if (i == 1) {
        # For i == 1, only the right bound cM is calculated as (rightBP) * r
        .y1$right_cM[i] <- (.y1[[3]][i] / 1e+6) * .y1[[4]][i]
        
        next
        
      }
      
      # For subsequent i, the right bound cM is calculated as (rightBP - leftBP) * r,
      # Left bound cM is previous right bound cM
      .y1$left_cM[i] <- .y1$right_cM[i - 1]
      .y1$right_cM[i] <- ( (.y1[[3]][i] - .y1[[2]][i]) / 1e+6) * .y1[[4]][i]
      
    }
    
    ## Cumulative sum both left_cM and right_cM
    .y1$left_cM <- cumsum(.y1$left_cM)
    .y1$right_cM <- cumsum(.y1$right_cM)
    
    
    ## Use .y1 for y values and right_BP for x values
    x_use <- c(0, .y1$right_BP)
    y_use <- c(0, .y1$right_cM)
    
    # x_out is the BP for which we want cM positions
    x_out <- .x[[2]]
    
    ## Interpolate x
    interp_cM <- approx(x = x_use, y = y_use, xout = x_out, rule = 2)
    
    # Return .x with cM positions
    cbind(.x, cM = interp_cM$y)
    
  }, SIMPLIFY = FALSE)
  
  # Bind rows and return
  x_interp <- do.call("rbind", interp_out)
  row.names(x_interp) <- row.names(x)
  
  return(x_interp)
  
}


#' Genotype probabilities given parent markers and the generation of selfing
#' 
#' @param f.gen The filial (F) generation (e.g. if \code{f.gen} is 3, this refers to an F_3 family).
#' @param r The recombination frequency. If passed, two-locus genotype probabilities are provided.
#' 
#' @return 
#' A \code{matrix} of conditional genotype probabilities for offspring marker genotypes (-1, 0, 1) (columns) given the marker genotypes of two parents (rows).
#' 
#' @export
#' 
#' 
geno_prob <- function(f.gen, r) {
  
  # Error checking
  stopifnot(is.numeric(f.gen))
  
  if (missing(r)) {
    # Matrix of probabilities
    probs <- rbind(
      both_het = c( (1 - 0.5^(1+f.gen)) / 2, 0.5^(1+f.gen), (1 - 0.5^(1+f.gen)) / 2 ),
      opp = c( (1 - 0.5^f.gen) / 2, 0.5^f.gen, (1 - 0.5^f.gen) / 2 ),
      two_one = c( 0.25 * (1 - 0.5^(f.gen-1)), 0.5^f.gen, 1 - ((0.25 * (1 - 0.5^(f.gen-1))) + (0.5^f.gen)) ),
      zero_one = rev(c( 0.25 * (1 - 0.5^(f.gen-1)), 0.5^f.gen, 1 - ((0.25 * (1 - 0.5^(f.gen-1))) + (0.5^f.gen)) )),
      zero_zero = c(1, 1e-50, 1e-50),
      two_two = c(1e-50, 1e-50, 1)
    )
    
    # Create a matrix
    colnames(probs) <- c(-1, 0, 1)
    row.names(probs) <- c("0|0", "-1|1", "0|1", "-1|0", "-1|-1", "1|1")
    
  } else {
    stopifnot(r >= 0, r <= 1)
    
    ## Two-locus diplotype probabilities
    ##
    ## From Broman2012'
    ##
    ## This assumes dimorphic parents at the two loci
    ## F1 = AA | BB or
    ## A || B
    ## A || B
    ##
    ## In other words, the parent genotypes are:
    ## 22 and 00, respectively
    ##
    ## Potential diplotypes are:
    ## 22
    ## 00
    ## 02
    ## 20
    ## 10
    ## 01
    ## 21
    ## 12
    ## 11
    ##
    ##
    
    # Set k to f.gen
    k <- f.gen
    
    ## Codify the probabilities from Broman2012
    probabilities <- c(
      "AA|AA" = ( 1 / (2 * (1+2*r)) ) - (( 0.5^(k+1) ) * ( 2 - ((1 - (2*r) + (2*r^2))^(k-1)) + ( (1 - 2*r)^k / (1 + 2*r) ) )),
      "AB|AB" = ( r / (1+2*r) ) - (( 0.5^(k+1) ) * (2 - ((1 - (2*r) + (2*r^2))^(k-1)) - ( ((1-2*r)^k) / (1+2*r) ))),
      "AA|AB" = ( 0.5^k ) * ( 1 - (1 - (2*r) + (2*r^2))^(k-1) ),
      "AA|BB" = ( 0.5^k ) * ( (1 - (2*r) + (2*r^2))^(k-1) + (1-2*r)^(k-1) ),
      "AB|BA" = ( 0.5^k ) * ( (1 - (2*r) + (2*r^2))^(k-1) - (1-2*r)^(k-1) )
    )
    
    ## Expand probabilities and return
    probs <- matrix(c(
      rep(probabilities[1], 2),
      rep(probabilities[2], 2),
      rep(probabilities[3], 4),
      probabilities[4] + probabilities[5]
    ), nrow = 1, dimnames = list("-1-1|11",  c("11", "-1-1", "-11", "1-1", "10", "01", "-10", "0-1", "00")))
    
  }
  
  return(probs)
  
}

