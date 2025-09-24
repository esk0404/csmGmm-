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

#------------------------------------------------------------------------------------
symm_fit_ind_single <- function(testStat, fixedMu, fixedPi) {
  # testStat: numeric vector of length K (e.g., c(Z1, Z2))
  # fixedMu: numeric vector of length K
  # fixedPi: vector of length 4 → prior over (0,0), (0,1), (1,0), (1,1)
  
  K <- length(testStat)
  
  # Define the four h^l configurations
  Hmat <- list(
    c(0, 0), # global null
    c(0, 1), # partial
    c(1, 0), # partial
    c(1, 1)  # full alt
  )
  
  # Compute density under each pattern
  dens <- sapply(1:4, function(idx) {
    hl_vec <- Hmat[[idx]]
    calc_dens_ind_hl(matrix(testStat, nrow = 1), fixedMu, pi = 0.5, hl_vec) # pi here controls sign mixture
  })
  
  # Weight by fixed priors
  weighted_dens <- dens * fixedPi
  
  # Total marginal probability of Z
  probZ <- sum(weighted_dens)
  
  # Null = configs with at least one zero
  null_indices <- c(1, 2, 3) # (0,0), (0,1), (1,0)
  probNull <- sum(weighted_dens[null_indices])
  
  # Local FDR
  lfdr <- probNull / probZ
  
  return(list(
    lfdr = lfdr,
    probZ = probZ,
    probNull = probNull
  ))
}


#------------------------------------------------------------------------------------
# Function to generate one test statistic (J=1)
generate_single_test <- function(mu, hl, pi = 0.5) {
  K <- length(mu)
  Z <- numeric(K)
  
  for (j in 1:K) {
    if (hl[j] == 0) {
      Z[j] <- rnorm(1, 0, 1)
    } else {
      # Sample sign according to probability pi for +mu and 1-pi for -mu
      sign_j <- sample(c(-1, 1), 1, prob = c(1 - pi, pi))
      Z[j] <- rnorm(1, mean = sign_j * mu[j], sd = 1)
    }
  }
  return(Z)
}

#----------------------------------------------------------------------------------------
fixedPi_90 <- c(0.03, 0.03, 0.04, 0.90)

run_simulation_single <- function(fixedMu, n_sim = 1000) {
  fdp_vec <- numeric(n_sim)
  power_vec <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    # Randomly choosing one h^l with probabiliy 0.03 for (0,0), 0.3 for (0,1), 0.04 for (1,0), and 0.9 for (1,1)
    state <- sample(1:4, 1, prob = fixedPi_90)
    
    # Map to hl vector
    hl_list <- list(
      c(0, 0), # global null
      c(0, 1), # partial
      c(1, 0), # partial
      c(1, 1)  # full alt
    )
    hl <- hl_list[[state]]
    
    # Generate one test statistic (J=1)
    testStat <- matrix(generate_single_test(fixedMu, hl), nrow = 1)
    true_label <- ifelse(any(hl != 0), 1, 0) # 1 = alt, 0 = null
    
    # Fit fixed model
    fit_result <- symm_fit_ind_single(testStat, fixedMu, fixedPi_90)
    lfdr <- fit_result$lfdr
    
    # Rejection rule: reject if lfdr <= 0.1
    reject <- (lfdr <= 0.1)
    
    if (reject) {
      if (true_label == 1) {
        fdp_vec[i] <- 0
        power_vec[i] <- 1
      } else {
        fdp_vec[i] <- 1
        power_vec[i] <- 0
      }
    } else {
      fdp_vec[i] <- 0
      power_vec[i] <- 0
    }
  }
  
  list(mean_fdp = mean(fdp_vec),
       power = mean(power_vec))
}

## Run function by varying effect sizes (mu)
mu_values <- seq(0, 8, by = 0.25)
results_single <- lapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)  # K = 2
  run_simulation_single(fixedMu)
})

mean_power_single <- sapply(results_single, function(res) res$power)
mean_fdp_single <- sapply(results_single, function(res) res$mean_fdp)

## Plot the power curve
plot(mu_values, mean_power_single, type = "l", pch = 19, col = "blue",
     ylim = c(0, 1), xlab = "Effect Size (mu)", ylab = "Power",
     main = "Power Curve (Single Testing, 90% Alt)")
abline(h = 0.8, col = "red", lty = 2)  # reference line for 80% power

## Plot the FDP curve
plot(mu_values, mean_fdp_single, type = "l", pch = 19, col = "blue",
     ylim = c(0, 1), xlab = "Effect Size (mu)", ylab = "FDP",
     main = "FDP Curve (Single Testing, 90% Alt)")
abline(h = 0.1, col = "red", lty = 2)  # reference FDR threshold

#----------------------------------------------------------------------------------
fixedPi_40 <- c(0.4, 0.1, 0.1, 0.4)

run_simulation_single <- function(fixedMu, n_sim = 1000) {
  fdp_vec <- numeric(n_sim)
  power_vec <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    # Randomly choosing one h^l with probability fixedPi_40
    state <- sample(1:4, 1, prob = fixedPi_40)
    
    # Map to hl vector
    hl_list <- list(
      c(0, 0), # global null
      c(0, 1), # partial
      c(1, 0), # partial
      c(1, 1)  # full alt
    )
    hl <- hl_list[[state]]
    
    # Generate one test statistic (J=1)
    testStat <- matrix(generate_single_test(fixedMu, hl), nrow = 1)
    true_label <- ifelse(any(hl != 0), 1, 0) # 1 = alt, 0 = null
    
    # Fit fixed model
    fit_result <- symm_fit_ind_single(testStat, fixedMu, fixedPi_90)
    lfdr <- fit_result$lfdr
    
    # Rejection rule: reject if lfdr <= 0.1
    reject <- (lfdr <= 0.1)
    
    if (reject) {
      if (true_label == 1) {
        fdp_vec[i] <- 0
        power_vec[i] <- 1
      } else {
        fdp_vec[i] <- 1
        power_vec[i] <- 0
      }
    } else {
      fdp_vec[i] <- 0
      power_vec[i] <- 0
    }
  }
  
  list(mean_fdp = mean(fdp_vec),
       power = mean(power_vec))
}

## Run function by varying effect sizes (mu)
mu_values <- seq(0, 8, by = 0.25)
results_single <- lapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)  # K = 2
  run_simulation_single(fixedMu)
})

mean_power_single <- sapply(results_single, function(res) res$power)
mean_fdp_single <- sapply(results_single, function(res) res$mean_fdp)

## Plot the power curve
plot(mu_values, mean_power_single, type = "l", pch = 19, col = "blue",
     ylim = c(0, 1), xlab = "Effect Size (mu)", ylab = "Power",
     main = "Power Curve (Single Testing, 40% Alt)")
abline(h = 0.8, col = "red", lty = 2)  # reference line for 80% power

## Plot the FDP curve
plot(mu_values, mean_fdp_single, type = "l", pch = 19, col = "blue",
     ylim = c(0, 1), xlab = "Effect Size (mu)", ylab = "FDP",
     main = "FDP Curve (Single Testing, 40% Alt)")
abline(h = 0.1, col = "red", lty = 2)  # reference FDR threshold

#--------------------------------------------------------------------------------
fixedPi_3 <- c(0.9, 0.04, 0.03, 0.03)

run_simulation_single <- function(fixedMu, n_sim = 1000) {
  fdp_vec <- numeric(n_sim)
  power_vec <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    # Randomly choosing one h^l with probability fixedPi_3
    state <- sample(1:4, 1, prob = fixedPi_3)
    
    # Map to hl vector
    hl_list <- list(
      c(0, 0), # global null
      c(0, 1), # partial
      c(1, 0), # partial
      c(1, 1)  # full alt
    )
    hl <- hl_list[[state]]
    
    # Generate one test statistic (J=1)
    testStat <- matrix(generate_single_test(fixedMu, hl), nrow = 1)
    true_label <- ifelse(any(hl != 0), 1, 0) # 1 = alt, 0 = null
    
    # Fit fixed model
    fit_result <- symm_fit_ind_single(testStat, fixedMu, fixedPi_90)
    lfdr <- fit_result$lfdr
    
    # Rejection rule: reject if lfdr <= 0.1
    reject <- (lfdr <= 0.1)
    
    if (reject) {
      if (true_label == 1) {
        fdp_vec[i] <- 0
        power_vec[i] <- 1
      } else {
        fdp_vec[i] <- 1
        power_vec[i] <- 0
      }
    } else {
      fdp_vec[i] <- 0
      power_vec[i] <- 0
    }
  }
  
  list(mean_fdp = mean(fdp_vec),
       power = mean(power_vec))
}

## Run function by varying effect sizes (mu)
mu_values <- seq(0, 8, by = 0.25)
results_single <- lapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)  # K = 2
  run_simulation_single(fixedMu)
})

mean_power_single <- sapply(results_single, function(res) res$power)
mean_fdp_single <- sapply(results_single, function(res) res$mean_fdp)

## Plot the power curve
plot(mu_values, mean_power_single, type = "l", pch = 19, col = "blue",
     ylim = c(0, 1), xlab = "Effect Size (mu)", ylab = "Power",
     main = "Power Curve (Single Testing, 3% Alt)")
abline(h = 0.8, col = "red", lty = 2)  # reference line for 80% power

## Plot the FDP curve
plot(mu_values, mean_fdp_single, type = "l", pch = 19, col = "blue",
     ylim = c(0, 1), xlab = "Effect Size (mu)", ylab = "FDP",
     main = "FDP Curve (Single Testing, 3% Alt)")
abline(h = 0.1, col = "red", lty = 2)  # reference FDR threshold

#---------------------------------------------------------------------------------
# Misspecify priors with fixedPi_3 as reference

run_simulation_single <- function(fixedMu, n_sim = 1000) {
  fdp_vec <- numeric(n_sim)
  power_vec <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    # Randomly choosing one h^l with probability fixedPi_3
    state <- sample(1:4, 1, prob = fixedPi_3)
    
    # Map to hl vector
    hl_list <- list(
      c(0, 0), # global null
      c(0, 1), # partial
      c(1, 0), # partial
      c(1, 1)  # full alt
    )
    hl <- hl_list[[state]]
    
    # Generate one test statistic (J=1)
    testStat <- matrix(generate_single_test(fixedMu, hl), nrow = 1)
    true_label <- ifelse(any(hl != 0), 1, 0) # 1 = alt, 0 = null
    
    # Fit fixed model
    fit_result <- symm_fit_ind_single(testStat, fixedMu, fixedPi_40)
    lfdr <- fit_result$lfdr
    
    # Rejection rule: reject if lfdr <= 0.1
    reject <- (lfdr <= 0.1)
    
    if (reject) {
      if (true_label == 1) {
        fdp_vec[i] <- 0
        power_vec[i] <- 1
      } else {
        fdp_vec[i] <- 1
        power_vec[i] <- 0
      }
    } else {
      fdp_vec[i] <- 0
      power_vec[i] <- 0
    }
  }
  
  list(mean_fdp = mean(fdp_vec),
       power = mean(power_vec))
}

## Run function by varying effect sizes (mu)
mu_values <- seq(0, 8, by = 0.25)
results_single <- lapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)  # K = 2
  run_simulation_single(fixedMu)
})

mean_power_single <- sapply(results_single, function(res) res$power)
mean_fdp_single <- sapply(results_single, function(res) res$mean_fdp)

## Plot the power curve
plot(mu_values, mean_power_single, type = "l", pch = 19, col = "blue",
     ylim = c(0, 1), xlab = "Effect Size (mu)", ylab = "Power",
     main = "Power Curve (Single Testing, 3% Alt)")
abline(h = 0.8, col = "red", lty = 2)  # reference line for 80% power

## Plot the FDP curve
plot(mu_values, mean_fdp_single, type = "l", pch = 19, col = "blue",
     ylim = c(0, 1), xlab = "Effect Size (mu)", ylab = "FDP",
     main = "FDP Curve (Single Testing, 3% Alt)")
abline(h = 0.1, col = "red", lty = 2)  # reference FDR threshold


