#' symm_fit_ind_fixed.R
#'
#' Fit the conditionally symmetric multidimensional Gaussian mixture model for sets of independent elements
#'
#' @param testStats J*K matrix of test statistics where J is the number of sets and K is number of elements in each set.
#'
#'
#' @return A list with the elements:
#' \item{fixedPi}{Vector of fixed pi values with length 3^K, holds the final mean parameters.}
#' \item{fixedMu}{Vector of fixed mu values with length K, holds final mixture proportions.}
#' \item{lfdrResults}{J*1 vector of all lfdr statistics.}
#' \item{Hmat}{Dataframe of all configurations including bl, sl, l, symAlt}
#' \item{probZ}{Vector of length J. For each row of testStats, the total probability (likelihood) under all configurations.}
#' \item{probNull}{Vector of length J. For each row of testStats, the total probability under the null configurations (configurations with at least one zero component in hl).}

#' @importFrom dplyr %>% mutate arrange filter select slice relocate
#' @import utils
#' @export
#' @examples
#' set.seed(0)
#' testStats <- rbind(null_stats, alt_1, alt_2, alt_3)
#' null_stats <- generate_alt_stats(900, fixedMu, pi, hl = c(0, 0))
#' alt_1 <- generate_alt_stats(40, fixedMu, pi, hl = c(1, 0)) 
#' alt_2 <- generate_alt_stats(30, fixedMu, pi, hl = c(0, 1)) 
#' alt_3 <- generate_alt_stats(30, fixedMu, pi, hl = c(1, 1)) 
#' fixedPi_3 <- c(0.9, 0.04, 0.03, 0.03)
#' fixedPi_90 <- c(0.03, 0.03, 0.04, 0.9)
#' fixedPi_40 <- c(0.4, 0.1, 0.1, 0.4)
#' fixedMu <- c(2, 2)
#' mu_values <- seq(0, 12, by = 0.25)
#' fixedMu <- rep(mu_val, 2)  : repeats scalar twice (e.g. when mu_values = 1, fixedMu = c(1,1))
#' results <- symm_fit_ind_fixed(testStats = testStats, fixedMu = fixedMu, fixedPi = fixedPi)


calc_dens_ind_hl <- function(Zmat, mu, pi, hl) {
  N <- nrow(Zmat)
  K <- ncol(Zmat)
  
  dens <- rep(1, N)
  
  for (j in 1:K) {
    dens_j <- if (hl[j] == 0) {
      dnorm(Zmat[,j], mean = 0, sd = 1)
    } else {
      pi * dnorm(Zmat[,j], mean = mu[j], sd = 1) +
        (1 - pi) * dnorm(Zmat[,j], mean = -mu[j], sd = 1)
    }
    
    dens <- dens * dens_j
  }
  
  return(dens)
}

library(dplyr)

symm_fit_ind_fixed <- function(testStats, fixedMu, fixedPi) {

  # number of composite null hypotheses
  J <- nrow(testStats)
  # number of dimensions
  K <- ncol(testStats)
  B <- 2^K - 1
  # number of hl configurations (excluding global null)
  L <- 3^K - 1

  # Only include the hl patterns that are actually used
  Hmat <- as.data.frame(matrix(c(
    0, 0,
    1, 0,
    0, 1,
    1, 1
  ), byrow = TRUE, ncol = K))

  # Precompute all densities
  densMat <- sapply(1:nrow(Hmat), function(idx) {
    # Extract h^l for this configuration
    hl_vec <- as.numeric(Hmat[idx, 1:K]) #Hmat[idx, ] gives one h^l configuration
    
    # For this problem, fixedPi is vector of 4 probabilities → grouped by pattern (regardless of signs)
    #fixedPi[idx] gives prior probability assigned to that particular h^l configuration
    calc_dens_ind_hl(testStats, fixedMu, fixedPi[idx], hl_vec)
  })
  
  # Enumerate all 9 signed configurations
  allHmat <- as.data.frame(expand.grid(c(-1, 0, 1), c(-1, 0, 1)))
  
  # Get the pattern index (row index of Hmat) for each signed config using abs()
  patternIndices <- apply(allHmat, 1, function(row) {
    match(paste0(abs(row), collapse=""), apply(Hmat, 1, function(x) paste0(x, collapse="")))
  })
  
  # Count how many signed configs belong to each pattern
  patternCounts <- table(patternIndices)
  
  # Reweight fixedPi across signed configs (each pattern's pi split equally)
  # Ex: if fixedPi = c(0.4, 0.1, 0.1, 0.4) then we weight h^l = (0,1), which is the third object in Hmat vector 
  #with prior 0.05 (the third object in fixedPi vector) because there are 2 h^ls with the same pattern (i.e. (0,1) and (0, -1)).
  expandedPi <- sapply(1:nrow(allHmat), function(i) {
    pIdx <- patternIndices[i]
    fixedPi[pIdx] / patternCounts[as.character(pIdx)]
  })
  
  # Expand densities using pattern index
  expandedDensMat <- sapply(1:nrow(allHmat), function(i) {
    densMat[, patternIndices[i]]
  })
  
  # Multiply by expanded priors
  weightedDensMat <- sweep(expandedDensMat, 2, expandedPi, `*`)
  
  # Compute probability of Z
  probZ <- rowSums(weightedDensMat)
  
  # Define null columns → any configuration with at least one 0
  nullCols <- which(apply(allHmat, 1, function(x) any(x == 0)))
  
  # Compute probNull
  probNull <- rowSums(weightedDensMat[, nullCols, drop=FALSE])
  
  # Local FDR
  lfdrResults <- probNull / probZ
  
  # Return info
  return(list(
    fixedPi = fixedPi,         # fixed pi values used
    fixedMu = fixedMu,         # fixed mu values used
    lfdrResults = lfdrResults, # computed local false discovery rates
    Hmat = Hmat,               # the configurations
    probZ = probZ,             # total prob
    probNull = probNull        # total null prob
  ))
}

#Create folder
#install.packages(c("usethis", "devtools"))
#usethis::create_package("C:/Users/yunsu/Downloads/Research/csmGmmFixed")



